/* 
 * File:   Mode.hpp
 * Author: matthewsupernaw
 *
 * Created on May 27, 2014, 1:29 PM
 */

#ifndef MODEL_HPP
#define	MODEL_HPP
#include "DataModule.hpp"

namespace ss {

    template<class REAL_T, class EVAL_T = REAL_T>
    class ModelBase {
        DataModule<REAL_T>* data;

    public:

        DataModule<REAL_T>* GetData() const {
            return data;
        }

        void SetData(DataModule<REAL_T>* data) {
            this->data = data;
        }




        virtual void Evaluate(EVAL_T &f) = 0;


    };
}




#endif	/* MODE_HPP */

