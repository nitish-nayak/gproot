# gproot

A set of simple ROOT-based classes for Gaussian Process based regression. These are mainly : 

* *GPPoint* : A container for points in 1D, 2D or 3D with a Mag() method describing the Euclidean distance between 2 points
* *GPKernel* : A base class that defines a kernel function over GPPoints along with any hyperparameters and is used as input to the gaussian process regression. Also calculates the kernel function derivatives. As such, only RBF and RationalQuadratic are currently implemented as they are the simplest kernel functions that are infinitely differentiable.
* *GPRegressor* : This handles the actual regression by taking in a list of GPPoints along with the corresponding data values. Needs a well defined kernel with some initial choice of hyperparameters and optionally a list of noise terms to add to the kernel diagonal. The hyperparameters are tuned by maximizing their log-marginal likelihood based on the fitted points and their data values. One can then obtain predictions at new points. 

Couple of things to note : 

* The hyperparameter tuning currently uses the BFGS algorithm for multi-dimensional optimization to mimic the popular GP implementation in scikit-learn. This is provided by the GSLMultiMin library through its interface in ROOT and thus requires a relatively modern version of ROOT (atleast 6+, I tested this on ROOT 6.18)
* The ROOT interface for the minimization is given by the base class IGradientFunctionMultiDim, which requires a well-defined gradient along with the function to be minimized. I don't know what this would mean for some classes of gaussian process kernels (for eg. Matern kernels) which don't require differentiability and in fact intrepret it as an additional hyperparameter to be tuned. Scikit-learn seems to handle it by not letting them float in the fit but instead, just minimizing for the rest of the parameters using only a small range of values for the differentiability parameter. Implementing this in ROOT requires some more thought. 

After compiling and adding it to the relevant include paths, as an illustration, one can do : 
```
std::vector<double> pars = {1., 1.};
RationalQuadraticKernel kern = RationalQuadraticKernel(pars);

GPRegressor reg = GPRegressor(kern);
reg.Fit(x, y);
reg.Predict(x_new);

TVectorD yhat = reg.PosteriorMean();
TVectorD ystd = reg.PosteriorStd();
```
to perform the regression and access their results within a ROOT macro. 

More details about Gaussian process and its implementation in scikit-learn can be found [here](https://scikit-learn.org/stable/modules/gaussian_process.html) and [here](https://scikit-learn.org/stable/auto_examples/gaussian_process/plot_gpr_noisy_targets.html). Information about different kinds of kernels can be found [here](https://www.cs.toronto.edu/~duvenaud/cookbook/). For Matern kernels specifically, see [here](https://scikit-learn.org/stable/modules/generated/sklearn.gaussian_process.kernels.Matern.html)
