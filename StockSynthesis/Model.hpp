/* 
 * File:   Mode.hpp
 * Author: matthewsupernaw
 *
 * Created on May 27, 2014, 1:29 PM
 */

#ifndef MODEL_HPP
#define	MODEL_HPP

#include <vector>

#include "DataModule.hpp"
#include "ModelFunctor.hpp"

namespace ss {

    template<class REAL_T, class EVAL_T = REAL_T>
    class ModelBase {
        typedef std::vector<ss::ModelFunctor<REAL_T, EVAL_T> > Functors;
        typedef typename std::vector<ss::ModelFunctor<REAL_T, EVAL_T> >::iterator FunctorsIterator;


        ss::DataModule<REAL_T>* data_m;
        Functors functors_m;
        EVAL_T value;

    public:

        operator EVAL_T() const {
            return this->Evaluate();
        }

        ss::DataModule<REAL_T>* GetData() const {
            return data_m;
        }

        void SetData(ss::DataModule<REAL_T>* data) {
            this->data_m = data;
        }

        void AddFunctor(ss::ModelFunctor<REAL_T, EVAL_T>* functor) {
            this->functors_m.push_back(functor);
        }

        virtual const EVAL_T Evaluate() {
            FunctorsIterator it;

            for (it = functors_m.begin(); it != functors_m.end(); ++it) {
                it->Evaluate();
            }

            return this->value;

        }

    private:



    };
}




#endif	/* MODE_HPP */

