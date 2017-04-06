# ccm
Convergent Cross Mapping

Implements Convergent Cross Mapping (CCM) (Sugihara et al. Science 2012)
2017/04/08 - Leonardo Barbosa

[x_r, y_r, x_o, y_o] = ccm(x,y,tau)

Parameters:

x, y : original time series
tau  : step between samples to construct the shadow manifold (default 1)
E    : size of the window to construct the manifold (default 10% original size)

Returns:

x_r, y_r : reconstructed time series
x_o, y_o : truncated original time series (so that first sample of x_o is
contemporaneous of x_r, i.e. reconstructed time series are missing N 
samples in the beggining, where N = 1+(E-1)*tau)

