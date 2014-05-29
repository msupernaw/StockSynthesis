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
        DataModule<REAL_T>* data_m;
        std::vector<ss::ModelFunctor<REAL_T, EVAL_T> > functors;

    public:

        DataModule<REAL_T>* GetData() const {
            return data_m;
        }

        void SetData(DataModule<REAL_T>* data) {
            this->data_m = data;
        }

        void AddFunctor(ss::ModelFunctor<REAL_T, EVAL_T>* functor) {

        }


        virtual void Evaluate(EVAL_T &f) = 0;


    };
}




#endif	/* MODE_HPP */

