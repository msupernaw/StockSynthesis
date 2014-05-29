/* 
 * File:   CatAtAge.hpp
 * Author: matthewsupernaw
 *
 * Created on May 29, 2014, 9:02 AM
 */

#ifndef CATCHATAGE_HPP
#define	CATCHATAGE_HPP

#include "../Model.hpp"
#include <vector>
#include <cmath>
#include <stdint.h>

namespace ss {

    template<class TT>
    const TT norm2(const std::vector<TT> &vect) {
        TT ret = TT(0.0);
        for (int i = 0; i < vect.size(); i++) {
            ret = ret + vect[i] * vect[i];
        }
        return ret;
    }

    template<class REAL_T>
    class CatchAtAgeData {
        int nyrs;
        int nages;
        std::vector<REAL_T> obs_catch_at_age;
        std::vector<REAL_T> effort;
        REAL_T M;
        std::vector<REAL_T> relwt;

    public:

        CatchAtAgeData() : nyrs(0), nages(0) {
        }

        REAL_T GetM() const {
            return M;
        }

        void SetM(REAL_T M) {
            this->M = M;
        }

        std::vector<REAL_T> GetEffort() const {
            return effort;
        }

        void SetEffort(std::vector<REAL_T> effort) {
            this->effort = effort;
        }

        int GetNages() const {
            return nages;
        }

        void SetNages(int nages) {
            this->nages = nages;
        }

        int GetNyrs() const {
            return nyrs;
        }

        void SetNyrs(int nyrs) {
            this->nyrs = nyrs;
        }

        std::vector<REAL_T> GetObs_catch_at_age() const {
            return obs_catch_at_age;
        }

        void SetObs_catch_at_age(std::vector<REAL_T> obs_catch_at_age) {
            this->obs_catch_at_age = obs_catch_at_age;
        }

        std::vector<REAL_T> GetRelwt() const {
            return relwt;
        }

        void SetRelwt(std::vector<REAL_T> relwt) {
            this->relwt = relwt;
        }

    };

    template<class REAL_T, class EVAL_T>
    class CatchAtAge;

    template<class REAL_T, class EVAL_T>
    class CatchAtAgeFunctor {
    public:
        virtual void Evaluate(CatchAtAge<REAL_T, EVAL_T> &model) = 0;


    };

    template<class REAL_T, class EVAL_T>
    class CatchAtAge {
        typedef std::vector<ss::CatchAtAgeFunctor<REAL_T, EVAL_T>* > Functors;
        typedef typename std::vector<ss::CatchAtAgeFunctor<REAL_T, EVAL_T>* >::iterator FunctorsIterator;

        uint32_t phase;
        uint32_t max_phase;
        CatchAtAgeData<REAL_T> data;

        ss::DataModule<REAL_T>* data_m;
        Functors functors_m;
        EVAL_T value;




        std::vector<EVAL_T> effort_devs;
        //Runtime


        std::vector<EVAL_T> F;
        std::vector<EVAL_T> Z;
        std::vector<EVAL_T> S;
        std::vector<EVAL_T> N;
        std::vector<EVAL_T> C;
        std::vector<EVAL_T> predicted_N;
        std::vector<EVAL_T> ratio_N;
        friend class CatchAtAgeFunctor<REAL_T, EVAL_T>;

        std::vector<EVAL_T*> active_parameters;
        std::vector<std::pair<EVAL_T*, uint32_t> >parameters;

        void AddParameter(EVAL_T* p, uint32_t phase) {
            parameters.push_back(std::pair<EVAL_T*, uint32_t > (p, phase));
        }

        void AddParameter(const std::vector<EVAL_T> &p, uint32_t phase) {
            for (int i = 0; i < p.size(); i++)
                parameters.push_back(std::pair<EVAL_T*, uint32_t > (&p[i], phase));
        }



    public:

        operator EVAL_T() const {
            return this->Evaluate();
        }

        REAL_T& GetValue() {
            return value;
        }

        void AddFunctor(ss::CatchAtAgeFunctor<REAL_T, EVAL_T>* functor) {
            this->functors_m.push_back(functor);
        }

        virtual const EVAL_T Evaluate() {
            FunctorsIterator it;

            for (int i = 0; i < this->functors_m.size(); i++) {
                this->functors_m[i]->Evaluate(*this);
            }

            return this->value;

        }

        const uint32_t Phase() {
            return this->phase;
        }

        void SetPhase(const uint32_t &phase) {
            this->phase = phase;
        }

        const uint32_t MaxPhase() {
            return this->max_phase;
        }

        void SetMaxPhase(const uint32_t &phase) {
            this->max_phase = phase;
        }

        CatchAtAgeData<REAL_T> GetData() const {
            return data;
        }

        void SetData(CatchAtAgeData<REAL_T> data) {
            this->data = data;
        }

        std::vector<EVAL_T> GetC() const {
            return C;
        }

        void SetC(std::vector<EVAL_T> C) {
            this->C = C;
        }

        std::vector<EVAL_T> GetF() const {
            return F;
        }

        void SetF(std::vector<EVAL_T> F) {
            this->F = F;
        }

        std::vector<EVAL_T> GetN() const {
            return N;
        }

        void SetN(std::vector<EVAL_T> N) {
            this->N = N;
        }

        std::vector<EVAL_T> GetS() const {
            return S;
        }

        void SetS(std::vector<EVAL_T> S) {
            this->S = S;
        }

        std::vector<EVAL_T> GetZ() const {
            return Z;
        }

        void SetZ(std::vector<EVAL_T> Z) {
            this->Z = Z;
        }

        std::vector<EVAL_T> GetEffort_devs() const {
            return effort_devs;
        }

        void SetEffort_devs(std::vector<EVAL_T> effort_devs) {
            this->effort_devs = effort_devs;
        }

        std::vector<EVAL_T> GetPredicted_N() const {
            return predicted_N;
        }

        void SetPredicted_N(std::vector<EVAL_T> predicted_N) {
            this->predicted_N = predicted_N;
        }

        std::vector<EVAL_T> GetRatio_N() const {
            return ratio_N;
        }

        void SetRatio_N(std::vector<EVAL_T> ratio_N) {
            this->ratio_N = ratio_N;
        }





    };

    template<class REAL_T, class EVAL_T>
    class MortalityAndSurvivability : public CatchAtAgeFunctor<REAL_T, EVAL_T> {
        EVAL_T log_q;
        std::vector<EVAL_T> log_sel;
        std::vector<EVAL_T> log_sel_coff;
    public:

        void Evaluate(CatchAtAge<REAL_T, EVAL_T> &model) {

            int i, j;
            // calculate the selectivity from the sel_coffs
            for (j = 0; j < model.GetData().GetNages() - 1; j++) {
                log_sel.at(j) = log_sel_coff.at(j);

            }

            if (model.GetData().GetNages() > 0) {
                // the selectivity is the same for the last two age classes
                log_sel.at(model.GetData().GetNages() - 1) = log_sel_coff.at(model.GetData().GetNages() - 2);
            }

            // This is the same as F(i,j)=exp(q)*effert(i)*exp(log_sel(j));
            //F = outer_prod(mfexp(log_q) * effort, mfexp(log_sel));
            for (i = 0; i < model.GetData().GetNyrs(); i++) {
                for (j = 0; j < model.GetData().GetNages(); j++) {
                    model.GetF().at(i * model.GetData().GetNages() + j)
                            = (std::exp(log_q) * model.GetEffort_devs().at(i)) * exp(log_sel.at(j));
                    //                std::cout<<F[i * nages + j]<<" "<<F[i * nages + j].wrt(this->log_q)<<"\n";
                    //                exit(0);
                }
            }

            if (model.Phase() == model.MaxPhase()) {//active(effort_devs)) {
                for (i = 0; i < model.GetData().GetNyrs(); i++) {
                    for (j = 0; j < model.GetData().GetNages(); j++) {
                        model.GetF().at(i * model.GetData().GetNages() + j)
                                = model.GetF().at(i * model.GetData().GetNages() + j) * std::exp(model.GetEffort_devs().at(i));
                    }
                }
            }
            //            
            //            
            //                        // get the total mortality
            for (i = 0; i < model.GetData().GetNyrs(); i++) {
                for (j = 0; j < model.GetData().GetNages(); j++) {
                    model.GetZ().at(i * model.GetData().GetNages() + j) =
                            model.GetF().at(i * model.GetData().GetNages() + j) + model.GetData().GetM();
                    model.GetS().at(i * model.GetData().GetNages() + j) =
                            std::exp(static_cast<REAL_T> (-1.0) * model.GetZ().at(i * model.GetData().GetNages() + j));
                }

            }
        }

    };

    template<class REAL_T, class EVAL_T>
    class NumbersAtAgeFunctor : public CatchAtAgeFunctor<REAL_T, EVAL_T> {
        EVAL_T log_popscale;
        std::vector<EVAL_T> log_relpop;
        std::vector<EVAL_T> log_initpop;

    public:

        void Evaluate(CatchAtAge<REAL_T, EVAL_T> &model) {
            int i, j;
            for (i = 0; i < log_initpop.size(); i++) {
                log_initpop.at(i) = log_relpop.at(i) + log_popscale;
            }

            for (i = 0; i < model.GetData().GetNyrs(); i++) {
                model.GetN().at(i * model.GetData().GetNages()) = std::exp(log_initpop.at(i));
            }

            for (j = 1; j < model.GetData().GetNages(); j++) {
                model.GetN().at(j) = std::exp(log_initpop.at((model.GetData().GetNyrs()) + j - 1));
            }
            //
            for (i = 0; i < model.GetData().GetNyrs() - 1; i++) {
                for (j = 0; j < model.GetData().GetNages() - 1; j++) {
                    model.GetN().at((i + 1) * model.GetData().GetNages() + (j + 1)) =
                            model.GetN().at(i * model.GetData().GetNages() + j) * model.GetS().at(i * model.GetData().GetNages() + j);
                }
            }

            model.GetValue() += (REAL_T) .01 * norm2<EVAL_T > (log_relpop);
        }
    };

    template<class REAL_T, class EVAL_T>
    class CatchFunctor : public CatchAtAgeFunctor<REAL_T, EVAL_T> {
    public:

        void Evaluate(CatchAtAge<REAL_T, EVAL_T> &model) {
            for (int i = 0; i < model.GetC().size(); i++) {
                model.GetC().at(i) = (model.GetF().at(i) / model.GetZ().at(i))*(((REAL_T) 1.0 - model.GetS().at(i)) * model.GetN().at(i));
            }
        }
    };

    template<class REAL_T, class EVAL_T>
    class CostFunctor : public CatchAtAgeFunctor<REAL_T, EVAL_T> {
    public:

        void Evaluate(CatchAtAge<REAL_T, EVAL_T> &model) {
    

            EVAL_T avg_F = (REAL_T) 0.0;
            for (int i = 0; i < model.GetF().size(); i++) {
                avg_F += model.GetF().at(i);
            }
            avg_F /= (REAL_T) model.GetF().size();

            if (model.Phase() == model.MaxPhase()) {

                // a very small penalty on the average fishing mortality
                model.GetValue() += (REAL_T) .001 * (std::log(avg_F / (REAL_T) .2) * std::log(avg_F / (REAL_T) .2));
            } else {
                model.GetValue() += (REAL_T) 1000. * (std::log(avg_F / (REAL_T) .2) * std::log(avg_F / (REAL_T) .2));
            }
            //
            EVAL_T sum;
            for (int i = 0; i < model.GetC().size(); i++) {
                sum += ((model.GetC().at(i) - model.GetData().GetObs_catch_at_age().at(i))*
                        (model.GetC().at(i) - model.GetData().GetObs_catch_at_age().at(i))) / ((REAL_T) 0.01 + model.GetC().at(i));
                //            std::cout << this->log_q << ":" << this->log_popscale << " " << f << "---" << sum.wrt(this->log_q) << "\n";
            }

            model.GetValue() += (REAL_T) 0.5 * REAL_T(model.GetC().size() + model.GetEffort_devs().size()) * std::log(sum + (REAL_T) 0.1 * norm2<EVAL_T > (model.GetEffort_devs()));

        }


    };



}




#endif	/* CATATAGE_HPP */

