/* 
 * File:   ModelFunctor.hpp
 * Author: matthewsupernaw
 *
 * Created on May 29, 2014, 7:58 AM
 */

#ifndef MODELFUNCTOR_HPP
#define	MODELFUNCTOR_HPP


//forward declaration for ModelBase.
template<class REAL_T, class EVAL_T>
class ss::ModelBase;




namespace ss {

    template<class REAL_T, class EVAL_T = REAL_T>
    class ModelFunctor {
        int rank_m;
        ss::ModelBase<REAL_T, EVAL_T>* model_m;
    public:

        ModelFunctor(ss::ModelBase<REAL_T, EVAL_T>* model) {

        }

        ss::ModelBase<REAL_T, EVAL_T>* GetModel() const {
            return model_m;
        }

        void SetModel(ss::ModelBase<REAL_T, EVAL_T>* model) {
            this->model_m = model;
        }

        int GetRank() const {
            return rank_m;
        }

        void SetRank(int rank) {
            this->rank_m = rank;
        }





    };

}


#endif	/* MODELFUNCTOR_HPP */

