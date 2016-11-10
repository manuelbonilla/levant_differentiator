/*!
 *	\file 	levant_differentiator.hpp
 *	\author	Manuel Bonilla (josemanuelbonilla@gmail.com)
 *	\date	November, 2016
 *	\brief This file contains an implementation of the exact differentiator \cite levant
*/

#ifndef LEVANT_DIFFERENTIATOR_HPP
#define LEVANT_DIFFERENTIATOR_HPP

#include <iostream>
#include <cmath>
#include <boost/numeric/odeint.hpp>

template <typename T>
class Levant {

	T L_; ///< Lipschitz constant L > 0
	int order_; ///< Order of the exact differentiator
	std::vector<T> lambda_; ///< Giains of the exact differentiator
	std::vector<T> x_; ///< output of the exact differentiator
	std::vector<T> x_filtered_; ///<Variable used in case that the filters is used
	std::vector<T> x_filtered_last_;
	boost::numeric::odeint::runge_kutta4< std::vector<T> > stepper; ///< integrator for the exact differentiator
	bool use_filter_; ///< Enable filter
	T alpha_; /*  alpha gain for the filter.
	           *  \f[
	           *  \alpha = (2*\pi*\Deltat*f_c) / (2*\pi*\Deltat*f_c + 1)
	           *  \f]
	           */

public:

	Levant();
	/**
	 * \brief      Sets the initial conditions for the exact differentiator
	 *
	 * \param[in]  x   Initial conditions
	 */
	void setInitialConditions(std::vector<T> x) {x_ = x;};
	/**
	 * \brief      returns the estimated functions
	 * 			x[0] contains the estimated input function
	 * 			x[n] contains the $n$-th derivative. First derivative is filtered if filter is activated
	 */

	std::vector<T> getX();
	/**
	 * \brief      make an step for the exact differentiator
	 *
	 * \param[in]  f     function to differentiate
	 * \param[in]  dt    sampling time
	 */
	void step(T f, T dt);
	/**
	 * \brief      Sets the parameters.
	 *
	 * \param[in]  order   Order of the exact differentiator
	 * \param[in]  L       Lipschitz constant L > 0 of the $n$-th derivative
	 * \param[in]  lambda  Gains for the exact differentiators. In the original paper Levant suggest {1.1,1.5,2,3,5,8} for the case of a differentiator of 1$^{st}$ to 5$^{th}$ order
	 */
	void setParameters(int order = 2, T L = 1.0, std::vector<T> lambda = {1.5, 1.1});
	/**
	 * \brief      enable filter if needed. By default the filter is not used. When the filter is enabled the first derivative of f is smoothed
	 *
	 * \param[in]  use_filter  The use filter
	 * \param[in]  f_c         Cut-off frequency for the Low-pass filter
	 * \param[in]  dt          sampling rate
	 */
	void setFilter(bool use_filter = false, T f_c = 1000, T dt = .01);
	/**
	 * \brief      this function explicitly implements the exact differentiator. The filter is a low-pass filter
	 *
	 * \param[in]  f     Integral of the exact differentiator
	 * \param [out]      States of the exact differentiator
	 * \param[in]  x_    Function to differentiate
	 */
	void Levant_function( const std::vector<T>  &f , std::vector<T>  &fp , const T x_ );
	/**
	 * This operator is useful when an external integrator is desired to use
	 */
	void operator() ( const std::vector<T>  &f , std::vector<T>  &fp , const T x_ )
	{
		Levant_function(f, fp, x_);
	}
};

/*! \class Levant
* \brief Implementation of the exact differentiator proposed by Arie Levant
*/

template <typename T>
Levant<T>::Levant()
{
	order_ = 1;
	L_ = 100.0;
	lambda_ = { 1.5, 1.1};
	x_ = std::vector<T>(order_ + 1, 0);
	x_filtered_ = std::vector<T>(order_ + 1, 0);
	x_filtered_last_ = std::vector<T>(order_ + 1, 0);
	use_filter_ = false;
}

template <typename T>
void Levant<T>::setParameters(int order, T L, std::vector<T> lambda)
{
	order_ = order;
	L_ = L;
	lambda_ = lambda;
	x_ = std::vector<T>(order_ + 1, 0);
	x_filtered_ = std::vector<T>(order_ + 1, 0);
	x_filtered_last_ = std::vector<T>(order_ + 1, 0);
	use_filter_ = false;
}

template < typename T>
void Levant<T>::step(T f, T dt)
{
	stepper.do_step(std::bind(&Levant<T>::Levant_function, &(*this), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3), this->x_, f, dt);
}

template < typename T>
void Levant<T>::setFilter(bool use_filter, T f_c, T dt)
{
	use_filter_ = use_filter;
	if (use_filter_)
	{
		alpha_ = (2 * M_PI * dt * f_c) / (2 * M_PI * dt * f_c + 1);
#if DEBUG
		std::cout << "Using filter Low-pas Filter with a cut-off frequency f_c = " << f_c << std::endl;
		std::cout << "alpha = " << alpha_ << std::endl;
#endif
	}
}

template < typename T>
std::vector<T> Levant<T>::getX()
{
	if (!use_filter_)
		return x_;

	for (int i = 0; i < order_; ++i)
	{
		if (i == 1)
		{
			x_filtered_[i] = alpha_ * x_[i] + (1.0 - alpha_) * x_filtered_last_[i];
		}
		else
		{
			x_filtered_[i]  = x_[i];
		}
		x_filtered_last_[i] = x_filtered_[i];
	}
#if DEBUG
	std::cout << "filtering. x = " << x_[1] << "\tx_filtered =" << x_filtered_[1] << std::endl;
#endif

	return x_filtered_;
}

template < typename T>
void Levant<T>::Levant_function( const std::vector<T>  &f , std::vector<T>  &fp , const T x_ )
{

	T order_T = (T)order_;
	fp[0] = -lambda_[0] * std::pow(L_, 1.0 / (order_T + 1.0))
	        * std::pow(std::abs(f[0] - x_), (order_T) / (order_T + 1))
	        * std::copysign(1.0, f[0] - x_)
	        + f[1];

	for (int i = 1; i < order_; ++i)
	{
		fp[i] = -lambda_[i] * std::pow(L_, 1.0 / (order_T - i + 1))
		        * std::pow(std::abs(f[i] - fp[i - 1]), (order_T - i) / (order_T - i  + 1))
		        * std::copysign(1.0, f[i] - fp[i - 1])
		        + f[i + 1];
	}

	fp[order_] = -lambda_[order_] * L_
	             * std::copysign(1.0, f[order_] - fp[order_ - 1]);
}
#endif //LEVANT_DIFFERENTIATOR_HPP