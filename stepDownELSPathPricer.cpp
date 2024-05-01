

#ifndef ANALYTICSOLVER_ANALYTICSOLVER_METHODS_MONTECARLO_STEPDOWNELSPATHPRICER_HPP_
#define ANALYTICSOLVER_ANALYTICSOLVER_METHODS_MONTECARLO_STEPDOWNELSPATHPRICER_HPP_

#include <utility>
#include <analyticsolver/methods/montecarlo/elspathpricer.hpp>
#include <analyticsolver/instruments/stepdownels.hpp>

namespace analyticsolver
{
	class StepDownELSPathPricer : public ELSPathPricer<MultiPath>
	{
	 public:
		StepDownELSPathPricer(StepDownELS::arguments args,
			ext::shared_ptr<YieldTermStructure> termStructure)
			: ELSPathPricer<MultiPath>(args),
			  termStructure_(std::move(termStructure))
		{
			Size nRedemption = args_.basicInfo.exerciseDates_.size();
			autocallProb_.resize(nRedemption);
			if (args_.couponInfo.hasCouponInfo_)
			{
				couponPayment_.resize(args_.couponInfo.couponDates_.size(),
					std::vector<Real>(args_.couponInfo.couponDates_[0].size(), 0.0));
			}
		}

		Real operator()(const MultiPath& path) const override;

	 private:
		ext::shared_ptr<YieldTermStructure> termStructure_;
	};

	inline Real StepDownELSPathPricer::operator()(const MultiPath& path) const
	{
		double payoff = 0.0;
		bool kiFlag = false;

		if (args_.redemptionInfo.isAutocalled_)
			return 0;

		if (args_.kiBarrierInfo.isHitKiBarrier_)
			kiFlag = true;

		///choose asset method
//		Path minPathByRate = makeSinglePath(path, Min);

		/// check lizard barrier
		std::vector<bool> lizardFlags(args_.basicInfo.exerciseDates_.size(), false);
		if (args_.lizardInfo.hasLizardInfo_)
		{
			for (Size i = 0; i < args_.basicInfo.exerciseDates_.size(); i++)
			{
				if (args_.basicInfo.exerciseDates_[i] > args_.basicInfo.valuationDate_
					&& args_.lizardInfo.checkLizardBarrier_[i] && !args_.lizardInfo.isHitLizardBarrier_[i])
				{
					Size interval = args_.basicInfo.exerciseDates_[i] - args_.basicInfo.valuationDate_;
					std::vector<Real> minValues;
					double minCheck;
					for (Size j = 0; j < args_.basicInfo.underlyingBasePrices_.size(); j++)
					{
						minValues.push_back(*std::min_element(path[j].begin(), path[j].begin() + interval));
					}
					minCheck = *std::min_element(minValues.begin(), minValues.end());
					if (minCheck > args_.lizardInfo.lizardBarrierLevel_[i])
					{
						lizardFlags[i] = true;
					}
				}
			}
		}

		///check KI barrier
		if (args_.kiBarrierInfo.hasKiBarrier_ && !kiFlag)
		{
			std::vector<Real> kiMinValues;
			double kiMinCheck;
			for (Size j = 0; j < args_.basicInfo.underlyingBasePrices_.size(); j++)
			{
				kiMinValues.push_back(*std::min_element(path[j].begin(), path[j].end()));
			}
			kiMinCheck = *std::min_element(kiMinValues.begin(), kiMinValues.end());
			if (kiMinCheck < args_.kiBarrierInfo.barrierLevel_)
				kiFlag = true;
		}

		/// 조기상환 / 만기상환 체크
		std::vector<Real> couponCheckValues;
		std::vector<Real> checkValues;
		double checkValue;
		double couponCheckValue;

		for (Size rdmIdx = 0; rdmIdx < args_.basicInfo.exerciseDates_.size(); rdmIdx++)
		{
			if (args_.basicInfo.exerciseDates_[rdmIdx] > args_.basicInfo.valuationDate_)
			{
				///쿠폰 있는경우**************************************************************************************
				if (args_.couponInfo.hasCouponInfo_)
				{
					Size couponPerRedemption = args_.couponInfo.couponDates_[rdmIdx].size();
					for (Size cIdx = 0; cIdx < couponPerRedemption; cIdx++)
					{
						if (args_.couponInfo.couponDates_[rdmIdx][cIdx]
							> args_.basicInfo.valuationDate_)///이거 계속 체크할 필요가 없음
						{
							///TODO daysMean 적용 방안 리서치
//							couponCheckValue = minPathByRate.at(args_.couponInfo.couponDates_[rdmIdx][cIdx]);
							couponCheckValues.clear();
							for (Size j = 0; j < args_.basicInfo.underlyingBasePrices_.size(); j++)
							{
								couponCheckValues.push_back(path[j].at(args_.couponInfo.couponDates_[rdmIdx][cIdx]));
							}
							couponCheckValue = *std::min_element(couponCheckValues.begin(), couponCheckValues.end());
							if (couponCheckValue > args_.couponInfo.strike_)
							{
								double df = termStructure_->discount(args_.couponInfo.couponDates_[rdmIdx][cIdx]);
								payoff += args_.couponInfo.couponRate_ * args_.basicInfo.notional_ * df;
								couponPayment_[rdmIdx][cIdx] += args_.couponInfo.couponRate_ *
									args_.basicInfo.notional_ * df;
							}
						}
					}
				}
				///쿠폰 있는경우**************************************************************************************
				checkValues.clear();
				for (Size j = 0; j < args_.basicInfo.underlyingBasePrices_.size(); j++)
				{
					checkValues.push_back(path[j].at(args_.basicInfo.exerciseDates_[rdmIdx]));
				}
				checkValue = *std::min_element(checkValues.begin(), checkValues.end());
//				checkValue = minPathByRate.at(args_.basicInfo.exerciseDates_[rdmIdx]);

				if (checkValue > args_.redemptionInfo.strikes_[rdmIdx])
				{
					double df = termStructure_->discount(args_.basicInfo.exerciseDates_[rdmIdx]);
					payoff += (1 + args_.redemptionInfo.couponRates_[rdmIdx]) *
						args_.basicInfo.notional_ * df;
					autocallProb_[rdmIdx] += 1;
					break;
				}
				else if (lizardFlags[rdmIdx])
				{
					double df = termStructure_->discount(args_.basicInfo.exerciseDates_[rdmIdx]);
					payoff = (1 + args_.lizardInfo.lizardDummyRate_[rdmIdx]) * args_.basicInfo.notional_ * df;
					autocallProb_[rdmIdx] += 1;
					break;
				}
				else if (rdmIdx == args_.basicInfo.exerciseDates_.size() - 1)
				{
					autocallProb_[rdmIdx] += 1;
					double df = termStructure_->discount(args_.basicInfo.exerciseDates_[rdmIdx]);
					if (!kiFlag)
					{
						if (args_.kiBarrierInfo.hasKiBarrier_)
						{
							// 낙인배리어가 있고 조기상환 안되고 낙인도 안됨
//							if (args_.maturityInfo.hasMatIrregularPayoff_)
//							{
//								/// 만기에 다른 손익구조가 적용되는 경우(ex. booster ELS)
//								payoff += matIrregularPayoff(checkValue) * args_.basicInfo.notional_ * df;
//								break;
//							}
							payoff += (1 + args_.redemptionInfo.couponRates_[rdmIdx]) *
								args_.basicInfo.notional_ * df;
							break;
						}
						else
						{
							//낙인 배리어가 없고 조기상환 안됨
//							if (args_.maturityInfo.hasMatIrregularPayoff_)
//							{
//								/// 만기에 다른 손익구조가 적용되는 경우(ex. booster ELS)
//								payoff += matIrregularPayoff(checkValue) * args_.basicInfo.notional_ * df;
//								break;
//							}
							payoff += (checkValue) * args_.basicInfo.notional_ * df;
							break;
						}
					}
					payoff += ((checkValue) * args_.kiBarrierInfo.kiInRate_ +
						args_.kiBarrierInfo.kidummyRate_) * args_.basicInfo.notional_ * df;

				}
			}
		}

		return payoff;
	}
}
#endif //ANALYTICSOLVER_ANALYTICSOLVER_METHODS_MONTECARLO_STEPDOWNELSPATHPRICER_HPP_


