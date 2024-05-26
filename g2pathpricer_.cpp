//
// Andrea_2023_05_11
//
#include <analyticsolver/methods/montecarlo/g2pathpricer.hpp>
#include <analyticsolver/cashflows/fixedratecoupon.hpp>
//extern std::chrono::system_clock::time_point start;
namespace analyticsolver
{
    G2PathPricer::G2PathPricer(TimeGrid timeGrid,
                               StructuredSwap::arguments args,
                               ext::shared_ptr<G2Model> model,
                               Size simNum,
                               Size polynomOrder,
                               LsmBasisSystem::PolynomType polynomType)
            : args_(std::move(args)), model_(std::move(model)), timeGrid_(std::move(timeGrid)), simIdx_(0),
              fixedValue_(0.0), polynomOrder_(polynomOrder), polynomType_(polynomType), coeffIdx_(0)
    {
        isCallable_ = !(args_.exercise == nullptr);
        deferredIdx_ = static_cast<Size>(timeGrid_.isDeferred());
        for (Size legIndex = 0; legIndex < args_.legs.size(); ++legIndex)
        {
            for (const auto& i: args_.legs[legIndex])
            {
                CashFlow& cf = *i;
                if (!cf.hasOccurred(args_.asOfDate))
                    break;
                exIdx_[legIndex]++;
            }
        }
        periodIdx_[0] = args_.fixedCoupons.size() - exIdx_[0];
        periodIdx_[1] = args_.structuredSchedule.size() - 1 - exIdx_[1];

        structuredLegValue_.resize(simNum, std::vector<Real>(periodIdx_[1], 0.0));
        periodDF_.resize(simNum, std::vector<DiscountFactor>(periodIdx_[1], 1.0));
        referenceRate_.resize(simNum);
        rateCoefficient_.reserve(simNum);
        setFixedLeg();
        if (isCallable_)
        {
            setLSMC();
        }
    }

    G2PathPricer::G2PathPricer(TimeGrid timeGrid,
                               StructuredSwap::arguments args,
                               const ext::shared_ptr<HullWhiteProcessArray>& processes,
                               Size simNum,
                               Size polynomOrder,
                               LsmBasisSystem::PolynomType polynomType)
            : args_(std::move(args)), timeGrid_(std::move(timeGrid)), simIdx_(0), fixedValue_(0.0),
              polynomOrder_(polynomOrder), polynomType_(polynomType), coeffIdx_(0)
    {
        isCallable_ = !(args_.exercise == nullptr);
        models_.resize(processes->size());
        for (Size i = 0; i < processes->size(); ++i)
            models_[i] = processes->process(i)->model();
        deferredIdx_ = static_cast<Size>(timeGrid_.isDeferred());
        for (Size legIndex = 0; legIndex < args_.legs.size(); ++legIndex)
        {
            for (const auto& i: args_.legs[legIndex])
            {
                CashFlow& cf = *i;
                if (!cf.hasOccurred(args_.asOfDate))
                    break;
                exIdx_[legIndex]++;
            }
        }
        periodIdx_[0] = args_.fixedCoupons.size() - exIdx_[0];
        periodIdx_[1] = args_.structuredSchedule.size() - 1 - exIdx_[1];

        structuredLegValue_.resize(simNum, std::vector<Real>(periodIdx_[1], 0.0));
        periodDF_.resize(simNum, std::vector<DiscountFactor>(periodIdx_[1], 1.0));
        referenceRate_.resize(simNum);
        rateCoefficient_.reserve(simNum);
        setFixedLeg();
        if (isCallable_)
        {
            setLSMC();
        }
    }

    void G2PathPricer::setLSMC()
    {
        v_ = LsmBasisSystem::pathBasisSystem(polynomOrder_, polynomType_);
        fixedLegIdx_.reserve(periodIdx_[0] - 1);
        for (Size i = 0; i < periodIdx_[0] - 1; ++i)
        {
            Date beforeEndDate;
            if( i == 0 )
            {
                beforeEndDate = args_.asOfDate;
            }
            else
            {
                beforeEndDate = args_.fixedEndDates[i + exIdx_[0] - 1];
            }

            for (Size j = 0; j < args_.exercise->dates().size(); ++j)
            {
                if (args_.exercise->date(j) <= args_.fixedEndDates[i + exIdx_[0]]
                    && args_.exercise->date(j) > beforeEndDate)
                {
                    fixedLegIdx_.push_back(static_cast<Integer>(i));
                    break;
                }
            }
        }
        exerciseIdx_.reserve(periodIdx_[1] - 1);
        for (Size i = 0; i < periodIdx_[1] - 1; ++i)
        {
            for (Size j = exIdx_[0]; j < args_.exercise->dates().size(); ++j)
            {
                if (args_.exercise->date(j) <= args_.structuredSchedule[i + 1 + exIdx_[1]]
                    && args_.exercise->date(j) > args_.structuredSchedule[i + exIdx_[1]])
                {
                    exerciseIdx_.push_back(i);
                    break;
                }
            }
        }
    }

    void G2PathPricer::setFixedLeg() const
    {
        if (isCallable_ && args_.exercise->zeroCallable())
        {
            fixedLegValue_.resize(periodIdx_[0], 0.0);
            cumulativeCoupon_.resize(args_.fixedCoupons.size(), 0.0);
            Rate fixedRate = ext::dynamic_pointer_cast<FixedRateCoupon>(args_.legs[0][0])->rate();
            const auto& fixedTimeFraction = args_.fixedTimeFraction;
            Real rateFactor = 1.0;
            Spread rateMargin = 0.0;
            cumulativeCoupon_[0] = args_.fixedCoupons[0];
            for (Size i = 0; i < args_.fixedCoupons.size(); ++i)
            {
                fixedLegValue_.back() += args_.exercise->compoundedValue(fixedRate,
                                                                         fixedTimeFraction.begin() + i, fixedTimeFraction.end(), rateFactor, rateMargin)
                                         * args_.nominal;

                if (i != 0)
                {
                    cumulativeCoupon_[i] = cumulativeCoupon_[i - 1]
                                           * args_.exercise->compoundFactor(fixedRate, rateFactor, rateMargin, fixedTimeFraction[i]);

                    cumulativeCoupon_[i] += args_.fixedCoupons[i];
                }
            }
        }
        else
        {
            fixedLegValue_.reserve(periodIdx_[0]);
            for (Size i = 0; i < args_.fixedCoupons.size(); ++i)
            {
                if (args_.asOfDate > args_.fixedPayDates[i])
                    continue;

                Real df = args_.discountCurve->discount(args_.fixedPayDates[i]);
                fixedValue_ += args_.fixedCoupons[i] * df;
                fixedLegValue_.push_back(args_.fixedCoupons[i]);
            }
        }
    }

    Real G2PathPricer::operator()(const MultiPath& multiPath) const
    {
        Size n = multiPath.pathSize();
        AS_REQUIRE(n > 0, "the path cannot be empty");
        Size assetNum = multiPath.assetNumber();
        AS_REQUIRE(assetNum > 0, "there must be some paths");

        const Date& valuationDate = args_.asOfDate;
        const std::vector<Size>& mandatoryIndexes = timeGrid_.mandatoryIndexes();
        referenceRate_[simIdx_].reserve(periodIdx_[0]);
        Size callIdx = 0;
        coeffIdx_ = 0;

        std::vector<Rate> refRates(2);
        Rate refRate, refRate2;
        std::vector<Rate> xPaths(2);
        DiscountFactor cumulativeDF = 1.0;
        Real value = 0.0;
        std::vector<StructuredCoupon::RateType> rateTypes;
        for (Size pIdx = 0; pIdx < periodIdx_[1]; ++pIdx)
        {
            auto cp = ext::dynamic_pointer_cast<StructuredCoupon>(args_.legs[1][pIdx + exIdx_[1]]);
            auto ibors = cp->indices();
            rateTypes = cp->rateType();

            Size deferredIdx = pIdx + deferredIdx_;
            xPaths[0] = multiPath[0][mandatoryIndexes[deferredIdx]];
            xPaths[1] = multiPath[1][mandatoryIndexes[deferredIdx]];
            const Date& startDate = args_.structuredResetDates[pIdx + exIdx_[1]];

            if (startDate < valuationDate)
            {
                refRate = args_.resetRate * cp->accrualPeriod();
            }
            else
            {
                if (models_.empty())
                    refRates = calculateRate(rateTypes, ibors, valuationDate, startDate, xPaths[0], xPaths[1]);   //
                else
                    refRates =   calculateRateUsingArray(rateTypes, ibors, valuationDate, startDate, xPaths);///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                refRate = refRates[0];
                refRate2 = refRates[1];
            }
            if (pIdx == 0 && exIdx_[1] != 0)
            {
                if (startDate < valuationDate && valuationDate < args_.structuredPayDates[exIdx_[1] - 1])
                {
                    Real prevCoupon = args_.resetRate;
                    prevCoupon *= args_.structuredAccrualTimes[exIdx_[1] - 1];
                    value += prevCoupon * args_.discountCurve->discount(args_.structuredPayDates[exIdx_[1] - 1]);
                    refRate = args_.nextResetRate * cp->accrualPeriod();
                }
            }

            /*
             *  coupon 변수에 스프레드를 저장한 후 콜옵션이 있을 경우 스프레드 정보를 넣어 주고
             *  그 다음 coupon 로직을 반영
             */
            Real coupon = refRate * cp->gearing()[0] + refRate2 * cp->gearing()[1];
            if (isCallable_)
            {
                const auto& callDates = args_.exercise->afterExerciseDates(args_.asOfDate);
                Size callSize = callDates.size();
                if (callIdx < callSize && callDates[callIdx] == cp->accrualStartDate())
                {
                    referenceRate_[simIdx_].push_back(coupon);
                    ++callIdx;
                }
            }
            coupon += cp->spread()[0] + cp->globalSpread();
            coupon *= cp->accrualPeriod();
            Real notionalExchange = (pIdx == periodIdx_[1] - 1) ? cp->nominalExchange() : 0.0;

            structuredLegValue_[simIdx_][pIdx] = args_.nominal * (notionalExchange + coupon);

            if (pIdx == 0 && deferredIdx_ != 0)
            {
                for (Size tIdx = 0; tIdx < mandatoryIndexes[deferredIdx_]; ++tIdx)
                    periodDF_[simIdx_][0] *= calculateDiscountFactor(timeGrid_[tIdx], timeGrid_[tIdx + 1],
                                                                     multiPath[0][tIdx], multiPath[1][tIdx], multiPath[0].randomNumber()[tIdx], tIdx);
            }
            for (Size tIdx = mandatoryIndexes[deferredIdx]; tIdx < mandatoryIndexes[deferredIdx + 1]; ++tIdx)
            {
                periodDF_[simIdx_][pIdx] *= calculateDiscountFactor(timeGrid_[tIdx], timeGrid_[tIdx + 1],
                                                                    multiPath[0][tIdx], multiPath[1][tIdx], multiPath[0].randomNumber()[tIdx], tIdx);
            }
            cumulativeDF *= periodDF_[simIdx_][pIdx];
            value += (notionalExchange + coupon) * cumulativeDF;
        }
        if (isCallable_ && args_.exercise->zeroCallable())
            fixedValue_ += fixedLegValue_.back() * cumulativeDF;
        simIdx_++;
        return args_.nominal * value;
    }

    void G2PathPricer::makeRateCoefficient(
            const std::vector<StructuredCoupon::RateType>& rateTypes,
            const std::vector<ext::shared_ptr<InterestRateIndex>>& ibors,
            const Date& valuationDate,
            const Date& startDate) const
    {
        std::vector<Rate> result(rateTypes.size(), 0.0);
        const DayCounter& dc = args_.structuredDayCount;
        for (Size i = 0; i < result.size(); ++i)
        {
            if (rateTypes[i] == StructuredCoupon::Spot)
            {
                const auto& iborIndex = ext::dynamic_pointer_cast<IborIndex>(ibors[i]);
                const Date& endDate = iborIndex->maturityDate(startDate);
                Time t0 = dc.yearFraction(valuationDate, startDate);
                Time T = dc.yearFraction(valuationDate, endDate);
                rateCoefficient_.push_back(model_->discountBondCoefficient(t0, T));
            }
            else if (rateTypes[i] == StructuredCoupon::Swap)
            {
                const auto& swapIndex = ext::dynamic_pointer_cast<SwapIndex>(ibors[i]);
                const Calendar& swapCalendar = swapIndex->fixingCalendar();
                BusinessDayConvention swapBDC = swapIndex->fixedLegConvention();
                const DayCounter& swapDC = swapIndex->dayCounter();
                const Period& swapCouponPeriod = swapIndex->fixedLegTenor();
                const Period& swapTenor = swapIndex->tenor();

                Time swapCouponTime = swapDC.yearFraction(valuationDate, startDate);
                g2rateCoefficient_.emplace_back(model_->g2discountBondCoefficient(0, swapCouponTime));

                Date prevCouponDate;
                Date couponEndDate = startDate;
                Period loop;
                while (loop < swapTenor)
                {
                    prevCouponDate = couponEndDate;
                    couponEndDate = swapCalendar.advance(couponEndDate, swapCouponPeriod, swapBDC);
                    Real couponTenor = swapDC.yearFraction(prevCouponDate, couponEndDate);
                    swapCouponTime = swapDC.yearFraction(valuationDate, couponEndDate);
//					rateCoefficient_.push_back(model_->discountBondCoefficient(0, swapCouponTime));
                    g2rateCoefficient_.emplace_back(model_->g2discountBondCoefficient(0, swapCouponTime));
                    loop += swapCouponPeriod;
                }
            }
        }
    }

    void G2PathPricer::makeRateCoefficient(
            const std::vector<StructuredCoupon::RateType>& rateTypes,
            const std::vector<ext::shared_ptr<InterestRateIndex>>& ibors,
            const Date& valuationDate,
            const Date& startDate,
            const std::vector<Rate>& xPaths) const
    {
        std::vector<Rate> result(rateTypes.size(), 0.0);
        const DayCounter& dc = args_.structuredDayCount;
        for (Size i = 0; i < result.size(); ++i)
        {
            if (rateTypes[i] == StructuredCoupon::Spot)
            {
                const auto& iborIndex = ext::dynamic_pointer_cast<IborIndex>(ibors[i]);
                const Date& endDate = iborIndex->maturityDate(startDate);
                Time t0 = dc.yearFraction(valuationDate, startDate);
                Time T = dc.yearFraction(valuationDate, endDate);
                rateCoefficient_.push_back(models_[i]->discountBondCoefficient(t0, T));
            }
            else if (rateTypes[i] == StructuredCoupon::Swap)
            {
                const auto& swapIndex = ext::dynamic_pointer_cast<SwapIndex>(ibors[i]);
                const Calendar& swapCalendar = swapIndex->fixingCalendar();
                BusinessDayConvention swapBDC = swapIndex->fixedLegConvention();
                const DayCounter& swapDC = swapIndex->dayCounter();
                const Period& swapCouponPeriod = swapIndex->fixedLegTenor();
                const Period& swapTenor = swapIndex->tenor();

                Date prevCouponDate;
                Date couponEndDate = startDate;
                Period loop;
                while (loop < swapTenor)
                {
                    prevCouponDate = couponEndDate;
                    couponEndDate = swapCalendar.advance(couponEndDate, swapCouponPeriod, swapBDC);
                    Real couponTenor = swapDC.yearFraction(prevCouponDate, couponEndDate);
                    Time swapCouponTime = swapDC.yearFraction(valuationDate, couponEndDate);
                    rateCoefficient_.push_back(models_[i]->discountBondCoefficient(0, swapCouponTime));
                    loop += swapCouponPeriod;
                }
            }
        }
    }

    std::vector<Rate> G2PathPricer::calculateRate(const std::vector<StructuredCoupon::RateType>& rateTypes,
                                                  const std::vector<ext::shared_ptr<InterestRateIndex>>& ibors,
                                                  const Date& valuationDate,
                                                  const Date& startDate,
                                                  Rate xPath,
                                                  Rate xPath2) const
    {
//        Rate xPath = xPaths[0];  //
//        Rate xPath2 = xPaths[1]; //
        if (simIdx_ == 0)
            makeRateCoefficient(rateTypes, ibors, valuationDate, startDate);
        std::vector<Rate> result(rateTypes.size(), 0.0);
        const DayCounter& dc = args_.structuredDayCount;
        for (Size i = 0; i < result.size(); ++i)
        {
            if (rateTypes[i] == StructuredCoupon::Spot)
            {
                const auto& iborIndex = ext::dynamic_pointer_cast<IborIndex>(ibors[i]);
                const Date& endDate = iborIndex->maturityDate(startDate);
                result[i] = (1 / model_->discountBond(rateCoefficient_[coeffIdx_],
                                                      dc.yearFraction(valuationDate, startDate), dc.yearFraction(valuationDate, endDate),
                                                      xPath, xPath2) - 1) / (dc.yearFraction(startDate, endDate));
                ++coeffIdx_;
            }
            else if (rateTypes[i] == StructuredCoupon::Swap)
            {
/*
                std::chrono::duration<double> end = std::chrono::system_clock::now() - start;
                if(simIdx_ == 8)
                    std::cout << "[Payoff_Layer] ===== CalculateRate Start: " << end.count() << std::endl;
*/
                const auto& swapIndex = ext::dynamic_pointer_cast<SwapIndex>(ibors[i]);
/*
                end = std::chrono::system_clock::now() - start;
                if(simIdx_ == 8)
                    std::cout << "[Payoff_Layer]swapIndex->maturityDate(startDate): " << end.count() << std::endl;
				const Date& endDate = swapIndex->maturityDate(startDate);
*/
/*
                end = std::chrono::system_clock::now() - start;
                if(simIdx_ == 8)
                    std::cout << "[Payoff_Layer]swapIndex->fixingCalendar(): " << end.count() << std::endl;
*/
                const Calendar& swapCalendar = swapIndex->fixingCalendar();
/*
                end = std::chrono::system_clock::now() - start;
                if(simIdx_ == 8)
                    std::cout << "[Payoff_Layer]swapIndex->fixedLegConvention(): " << end.count() << std::endl;
*/
                BusinessDayConvention swapBDC = swapIndex->fixedLegConvention();
/*				if(simIdx_ == 0)
                    cmsFixedBDC_.emplace_back(swapIndex->fixedLegConvention());
                BusinessDayConvention swapBDC = cmsFixedBDC_[fixedBDCIdx_];
                ++fixedBDCIdx_;
*/
/*
                end = std::chrono::system_clock::now() - start;
                if(simIdx_ == 8)
                    std::cout << "[Payoff_Layer]swapIndex->dayCounter(): " << end.count() << std::endl;
*/
                const DayCounter& swapDC = swapIndex->dayCounter();
/*
                end = std::chrono::system_clock::now() - start;
                if(simIdx_ == 8)
                    std::cout << "[Payoff_Layer]swapIndex->fixedLegTenor(): " << end.count() << std::endl;
*/
                const Period& swapCouponPeriod = swapIndex->fixedLegTenor();
/*
                end = std::chrono::system_clock::now() - start;
                if(simIdx_ == 8)
                    std::cout << "[Payoff_Layer]swapIndex->tenor(): " << end.count() << std::endl;
*/
                const Period& swapTenor = swapIndex->tenor();

                Date prevCouponDate;
                Date couponEndDate = startDate;
                Real zeroBondSum = 0.0;
/*
                end = std::chrono::system_clock::now() - start;
                if(simIdx_ == 8)
                    std::cout << "[Payoff_Layer]model_->g2discountBond(..., startDate, ...): " << end.count() << std::endl;
                Rate swapRate =
                    model_->discountBond(0, dc.yearFraction(valuationDate, startDate), xPath, xPath2)
                        - model_->discountBond(0, dc.yearFraction(valuationDate, endDate), xPath, xPath2);
*/
                Rate swapRate = model_->g2discountBond(g2rateCoefficient_[coeffIdx_], 0,
                                                       swapDC.yearFraction(valuationDate, startDate), xPath, xPath2);
                ++coeffIdx_;

/*
              end = std::chrono::system_clock::now() - start;
              if(simIdx_ == 8)
                  std::cout << "[Payoff_Layer]swapRate Loop: " << end.count() << std::endl;
*/
                Period loop;
                while (loop < swapTenor)
                {
                    prevCouponDate = couponEndDate;
                    couponEndDate = swapCalendar.advance(couponEndDate, swapCouponPeriod, swapBDC);
                    Real couponTenor = swapDC.yearFraction(prevCouponDate, couponEndDate);

                    zeroBondSum +=
                            model_->g2discountBond(g2rateCoefficient_[coeffIdx_], 0,
                                                   swapDC.yearFraction(valuationDate, couponEndDate), xPath, xPath2) * couponTenor;
//                        model_->discountBond(rateCoefficient_[coeffIdx_], 0,
                    loop += swapCouponPeriod;
                    ++coeffIdx_;
                }
/*
                end = std::chrono::system_clock::now() - start;
                if(simIdx_ == 8)
                    std::cout << "[Payoff_Layer]swapRate calc end: " << end.count() << std::endl;
//                std::cout << " =========== endDate, underlyingSwapMaturityDate: " << endDate << "," << couponEndDate << std::endl;
*/
//                Size tempCoffIdx = coeffIdx_-1;
                swapRate -= model_->g2discountBond(g2rateCoefficient_[coeffIdx_-1], 0,
                                                   swapDC.yearFraction(valuationDate, couponEndDate), xPath, xPath2);

/*
                Rate tempSwapRate = model_->discountBond(0, dc.yearFraction(valuationDate, startDate), xPath, xPath2)
                                    - model_->g2discountBond(g2rateCoefficient_[coeffIdx_-1], 0,
                                                      swapDC.yearFraction(valuationDate, couponEndDate), xPath, xPath2);

*/
/*
                end = std::chrono::system_clock::now() - start;
                if(simIdx_ == 8)
                    std::cout << "[Payoff_Layer]swapRate /= zeroBondSum: " << end.count() << std::endl;
*/
                swapRate /= zeroBondSum;
                result[i] = swapRate;
            }
        }
        return result;
    }

    std::vector<Rate> G2PathPricer::calculateRateUsingArray(const std::vector<StructuredCoupon::RateType>& rateTypes,
                                                            const std::vector<ext::shared_ptr<InterestRateIndex>>& ibors,
                                                            const Date& valuationDate,
                                                            const Date& startDate,
                                                            const std::vector<Rate>& xPaths) const
    {
        if (simIdx_ == 0)
            makeRateCoefficient(rateTypes, ibors, valuationDate, startDate, xPaths);
        std::vector<Rate> result(rateTypes.size(), 0.0);
        const DayCounter& dc = args_.structuredDayCount;
        for (Size i = 0; i < result.size(); ++i)
        {
            if (rateTypes[i] == StructuredCoupon::Spot)
            {
                const auto& iborIndex = ext::dynamic_pointer_cast<IborIndex>(ibors[i]);
                const Date& endDate = iborIndex->maturityDate(startDate);
                result[i] = (1 / models_[i]->discountBond(rateCoefficient_[coeffIdx_],
                                                          dc.yearFraction(valuationDate, startDate), dc.yearFraction(valuationDate, endDate), xPaths[i]) - 1)
                            / (dc.yearFraction(startDate, endDate));
                ++coeffIdx_;
//				result[i] = (1 / models_[i]->discountBond(dc.yearFraction(valuationDate, startDate),
//					dc.yearFraction(valuationDate, endDate), xPaths[i]) - 1) /
//					(dc.yearFraction(startDate, endDate));
            }
            else if (rateTypes[i] == StructuredCoupon::Swap)
            {
                const auto& swapIndex = ext::dynamic_pointer_cast<SwapIndex>(ibors[i]);
                const Date& endDate = swapIndex->maturityDate(startDate);

                const Calendar& swapCalendar = swapIndex->fixingCalendar();
                BusinessDayConvention swapBDC = swapIndex->fixedLegConvention();
                const DayCounter& swapDC = swapIndex->dayCounter();
                const Period& swapCouponPeriod = swapIndex->fixedLegTenor();
                const Period& swapTenor = swapIndex->tenor();

                Date prevCouponDate;
                Date couponEndDate = startDate;
                Real zeroBondSum = 0.0;
                Rate swapRate =
                        models_[i]->discountBond(0, dc.yearFraction(valuationDate, startDate), xPaths[i])
                        - models_[i]->discountBond(0, dc.yearFraction(valuationDate, endDate), xPaths[i]);

                Period loop;
                while (loop < swapTenor)
                {
                    prevCouponDate = couponEndDate;
                    couponEndDate = swapCalendar.advance(couponEndDate, swapCouponPeriod, swapBDC);
                    Real couponTenor = swapDC.yearFraction(prevCouponDate, couponEndDate);
                    zeroBondSum +=
                            models_[i]->discountBond(rateCoefficient_[coeffIdx_], 0,                                                                                           ///////////////////////////////////// 7th hit
                                                     swapDC.yearFraction(valuationDate, couponEndDate),xPaths[i]) * couponTenor;
//					zeroBondSum +=
//						models_[i]->discountBond(0, swapDC.yearFraction(valuationDate, couponEndDate),xPaths[i]) * couponTenor;
                    loop += swapCouponPeriod;
                    ++coeffIdx_;
                }
                swapRate /= zeroBondSum;
                result[i] = swapRate;
            }
        }
        return result;
    }

    DiscountFactor G2PathPricer::calculateDiscountFactor(Time t,
                                                         Time T,
                                                         Real randomNumber,
                                                         Rate xPath1,
                                                         Real xPath2,
                                                         Size tIdx) const
    {
        bool afterSimulation = (simIdx_ != 0);
        DiscountFactor result = (models_.empty()) ?
                                model_->discount(t, T, xPath1, xPath2, randomNumber, afterSimulation, tIdx) :
                                models_[0]->discount(t, T, xPath1, randomNumber, afterSimulation, tIdx);
        return result;
    }

    void G2PathPricer::setRefRate(const std::vector<Rate>& rates) const
    {
        const auto& callDates = args_.exercise->afterExerciseDates(args_.asOfDate);
        Size callSize = callDates.size();
        referenceRate_[simIdx_].resize(callDates.size(), 0.0);
        for (Size pIdx = 1, callIdx = 0; pIdx < periodIdx_[1]; ++pIdx)
        {
            AS_REQUIRE( pIdx + exIdx_[1] <= args_.structuredResetDates.size(), "Size of period is greater than size of reset date." )
            if (callIdx < callSize
                && callDates[callIdx] <= args_.structuredResetDates[pIdx + exIdx_[1]]
                && callDates[callIdx] > args_.structuredResetDates[pIdx + exIdx_[1] - 1])
            {
                referenceRate_[simIdx_][callIdx] = rates[pIdx];
                ++callIdx;
            }
        }
    }

    Real G2PathPricer::optionValue() const
    {
        AS_REQUIRE(!referenceRate_[0].empty(), "Reference Rates at Call Dates are Null");

        if (args_.asOfDate >= args_.exercise->dates().back())
            return 0.0;

        if (args_.exercise->zeroCallable())
            return zeroCallValue();

        boost::scoped_array<Array> coeff_(new Array[periodIdx_[1] - 2]);
        Array optionValue(simIdx_, 0.0), exercise(simIdx_, 0.0);
        std::vector<Real> x;
        std::vector<Real> y;

        std::vector<Integer> payer(2);
        Integer callType = static_cast<Integer>(args_.exercise->callType());
        payer[0] = args_.payer[0] * callType;
        payer[1] = args_.payer[1] * callType;

        for (Size sIdx = 0; sIdx < simIdx_; sIdx++)
        {
            Real lastPayment = -payer[0] * fixedLegValue_.back();
            for (Size pIdx = periodIdx_[1] - 1; pIdx > exerciseIdx_.back(); --pIdx)
            {
                lastPayment -= payer[1] * structuredLegValue_[sIdx][pIdx];
                lastPayment *= periodDF_[sIdx][pIdx];
            }
            exercise[sIdx] = lastPayment;
            if (lastPayment > 0.0)
                optionValue[sIdx] = lastPayment;
        }

//		Size fIdx = fixedLegValue_.size() - 1;
        Size rIdx = exerciseIdx_.size();
        auto fIter = fixedLegIdx_.rbegin();
        for (auto iter = exerciseIdx_.rbegin(); iter != exerciseIdx_.rend() - 1; ++iter)
        {
            x.clear();
            y.clear();
            --rIdx;
//			--fIdx;
            for (Size sIdx = 0; sIdx < simIdx_; ++sIdx)
            {
                /*DiscountFactor cumulativeDF = 1.0;
                for (Size iterIdx = *iter; iterIdx > *(iter + 1); --iterIdx)
                {
                    Real fixedValue = (iterIdx == *iter) ? fixedLegValue_[fIdx] : 0.0;
                    exercise[sIdx] -=
                        args_.payer[1] * structuredLegValue_[sIdx][iterIdx] + args_.payer[0] * fixedValue;
                    exercise[sIdx] *= periodDF_[sIdx][iterIdx];
                    cumulativeDF *= periodDF_[sIdx][iterIdx];
                }
                optionValue *= cumulativeDF;
                for (Size iterIdx = *iter; iterIdx > *(iter + 1); --iterIdx)
                {
                    Real fixedValue = (iterIdx == *iter) ? fixedLegValue_[fIdx] : 0.0;
                    exercise[sIdx] -=
                        args_.payer[1] * structuredLegValue_[sIdx][iterIdx] + args_.payer[0] * fixedValue;
                    exercise[sIdx] *= periodDF_[sIdx][iterIdx];
                    optionValue[sIdx] *= periodDF_[sIdx][iterIdx];
                }*/
                Integer fIdx = *fIter;
                for (Size iterIdx = *iter; iterIdx > *(iter + 1); --iterIdx)
                {
                    Real fixedValue = (fIdx <= *(fIter + 1)) ? 0.0 : fixedLegValue_[fIdx];
                    exercise[sIdx] -=
                            payer[1] * structuredLegValue_[sIdx][iterIdx] + payer[0] * fixedValue;
                    exercise[sIdx] *= periodDF_[sIdx][iterIdx];
                    optionValue[sIdx] *= periodDF_[sIdx][iterIdx];
                    --fIdx;
                }
                if (exercise[sIdx] > 0.0)
                {
                    x.push_back(referenceRate_[sIdx][rIdx - 1]);
                    y.push_back(optionValue[sIdx]);
                }
            }
            ++fIter;

            if (v_.size() <= x.size())
                coeff_[*iter - 1] = GeneralLinearLeastSquares(x, y, v_).coefficients();
            else
                coeff_[*iter - 1] = Array(v_.size(), 0.0);

            for (Size sIdx = 0, idx = 0; sIdx < simIdx_; ++sIdx)
            {
                if (exercise[sIdx] > 0.0)    // exerices value : 행사시 가치는 이전 payoff도 다 가져와야 함.
                {
                    Real continuationValue = 0.0;
                    for (Size vIdx = 0; vIdx < v_.size(); ++vIdx)
                        continuationValue += coeff_[*iter - 1][vIdx] * v_[vIdx](x[idx]);
                    if (continuationValue < exercise[sIdx])
                        optionValue[sIdx] = exercise[sIdx];
                    ++idx;
                }
            }
        }

        Real value = 0.0;
        for (Size sIdx = 0; sIdx < simIdx_; sIdx++)
        {
            for (Integer pIdx = exerciseIdx_.front(); pIdx >= 0; pIdx--)
            {
                value += optionValue[sIdx] * periodDF_[sIdx][pIdx];
            }
        }
        return value / simIdx_;
    }

    Real G2PathPricer::zeroCallValue() const
    {
        boost::scoped_array<Array> coeff_(new Array[periodIdx_[1] - 2]);
        Array optionValue(simIdx_, 0.0), exercise(simIdx_, 0.0);
        std::vector<Real> x;
        std::vector<Real> y;

//		Size cIdx = cumulativeCoupon_.size() - 2;    // 일부러 인덱스 하나 적게 설정
        for (Size sIdx = 0; sIdx < simIdx_; sIdx++)
        {
            Real lastPayment = -args_.payer[0] * cumulativeCoupon_.back();
            for (Size pIdx = periodIdx_[1] - 1; pIdx > exerciseIdx_.back(); --pIdx)
            {
                lastPayment -= args_.payer[1] * structuredLegValue_[sIdx][pIdx];
                lastPayment *= periodDF_[sIdx][pIdx];
            }
            lastPayment += args_.payer[0] * cumulativeCoupon_[exerciseIdx_.back() + exIdx_[0]];
            if (lastPayment > 0.0)
                optionValue[sIdx] = lastPayment;
        }

        Size rIdx = exerciseIdx_.size();
        auto fIter = fixedLegIdx_.rbegin() + 1;
        for (auto iter = exerciseIdx_.rbegin(); iter != exerciseIdx_.rend() - 1; ++iter)
        {
            x.clear();
            y.clear();
            --rIdx;
            for (Size sIdx = 0; sIdx < simIdx_; ++sIdx)
            {
                exercise[sIdx] = -args_.payer[0] * cumulativeCoupon_.back();
                DiscountFactor cumulativeDF = 1.0;
                Integer fIdx = *fIter;
                for (Size pIdx = periodIdx_[1] - 1; pIdx > *(iter + 1); --pIdx)
                {
                    exercise[sIdx] -= args_.payer[1] * structuredLegValue_[sIdx][pIdx];
                    exercise[sIdx] *= periodDF_[sIdx][pIdx];
                    if (pIdx <= *iter)
                        cumulativeDF *= periodDF_[sIdx][pIdx];
                }
                exercise[sIdx] += args_.payer[0] * cumulativeCoupon_[fIdx + exIdx_[0]];
                optionValue[sIdx] *= cumulativeDF;
                if (exercise[sIdx] > 0.0)
                {
                    x.push_back(referenceRate_[sIdx][rIdx - 1]);
                    y.push_back(optionValue[sIdx]);
                }
            }
            ++fIter;

            if (v_.size() <= x.size())
                coeff_[*iter - 1] = GeneralLinearLeastSquares(x, y, v_).coefficients();
            else
                coeff_[*iter - 1] = Array(v_.size(), 0.0);

            for (Size sIdx = 0, idx = 0; sIdx < simIdx_; ++sIdx)
            {
                if (exercise[sIdx] > 0.0)    // exerices value : 행사시 가치는 이전 payoff도 다 가져와야 함.
                {
                    Real continuationValue = 0.0;
                    for (Size vIdx = 0; vIdx < v_.size(); ++vIdx)
                        continuationValue += coeff_[*iter - 1][vIdx] * v_[vIdx](x[idx]);
                    if (continuationValue < exercise[sIdx])
                        optionValue[sIdx] = exercise[sIdx];
                    ++idx;
                }
            }
        }

        Real value = 0.0;
        for (Size sIdx = 0; sIdx < simIdx_; sIdx++)
        {
            for (Integer pIdx = exerciseIdx_.front(); pIdx >= 0; pIdx--)
            {
                value += optionValue[sIdx] * periodDF_[sIdx][pIdx];
            }
        }
        return value / simIdx_;
    }

    Disposable<Array> G2PathPricer::exercisedValue() const
    {
        AS_REQUIRE(args_.resetRate != Null<Rate>(), "reset rate dose not provided");
        Array value(2);
        auto cp = ext::dynamic_pointer_cast<StructuredCoupon>(args_.legs[1][exIdx_[1]]);
        Real coupon = args_.nominal * cp->accrualPeriod() * args_.resetRate;
        value[1] = coupon * args_.discountCurve->discount(args_.structuredPayDates[exIdx_[1]]);

        if (args_.exercise->zeroCallable())
        {
            value[0] = cumulativeCoupon_[exIdx_[0]];
        }
        else
        {
            value[0] = fixedLegValue_.front();
        }

        return value;
    }
}