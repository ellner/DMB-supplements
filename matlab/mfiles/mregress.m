function [Coefficients, S_err, XTXI, R_sq, F_val, Coef_stats, Y_hat, residuals, covariance] ...
          = mregress(Y, X, INTCPT)

% MREGRESS  Performs multiple linear regression analysis of X (independent) on Y (dependent).
%
%           Usage:
%
%           [Coefficients, S_err, XTXI, R_sq, F_val, Coef_stats, Y_hat, residuals, covariance] ... 
%                                                     = mregress(Y, X, INTCPT)
%
%	    INTCPT = 1; include a y-intercept in the model
%           INTCPT = 0; DO NOT include a y-intercept in the model
%
%           Returns:
%
%           Coefficients         - Regression coefficients
%           S_err                - Standard error of estimate
%           XTXI                 - inverse of X' * X
%           R_sq                 - R-squared
%           F_val                - F-value for the regression and siginificance level (p-value for F)
%           Coef_stats           - Coefficients with their standard deviations, T-values, and p-values
%           Y_hat                - Fitted values
%           residuals            - Residuals
%	    covariance           - Covariance matrix ( XTXI * S_err^2 )
%
%
% G. Anthony Reina
% Motor Control Lab
% The Neurosciences Institute
% Created: 4 Aug 1998
% Last Update: 10/8/1998 by GAR
%
% Please note that for the case when the intercept of the model equals zero, the
% definition of R-squared and the F-statistic change mathematically. For a linear
% model containing the y-intercept, R-squared refers to the amount of variance around
% the mean of the dependent variable (y) which is explained by the variance around the
% mean of the independent variables (x). For a linear model NOT containing the
% y-intercept, R-squared measures the amount of variance around ZERO of the dependent
% variable which is explained by the variance around ZERO of the independent variable.
% If the same equation for R-squared is used for both with and without a y-intercept
% (namely R-squared = [Sum of Squares of the Regression] / [ Total sum of the squares]),
% then R-squared may be a NEGATIVE value for some data. For this reason, 
% this subroutine will calculate R-squared using the total un-corrected sum of the
% squares. In effect, this approach avoids negative R-squares but may lack any
% meaningful interpretation for the "goodness-of-fit" in the model. It has been
% suggested by some texts that a more useful approach is to always use the case
% where y-intercept is included in the model. However, as with all statistical
% analyses, it is prudent to simply be aware of what your data "looks" like and
% what your statistical tools are actually measuring in order to generate a useful
% analysis.
%
% For further reading on regression through the origin (i.e. without a y-intercept),
% please refer to:
%
%    Neter J, Kutner MH, Nachtsheim CJ, and Wasserman W. "Applied Linear 
%             Statistical Models" 4th ed. Irwin publishing (Boston, 1996), pp 159-163.
%
%    Myers R, "Classical and Modern Regression with Applications" Duxbury Press
%             (Boston, 1986), p. 30.
%
%

if (nargin < 2)
   error('mregress requires at least 2 input variables. Type ''help mregress''.');
end

if (nargin == 2) 
   INTCPT = 0;
end

% Check that independent (X) and dependent (Y) data have compatible dimensions
% ----------------------------------------------------------------------------
[n_x, k] = size(X);
[n_y,columns] = size(Y);

if n_x ~= n_y, 
    error('The number of rows in Y must equal the number of rows in X.'); 
end 

if columns ~= 1, 
    error('Y must be a vector, not a matrix'); 
end

n = n_x;

%  Solve for the regression coefficients using ordinary least-squares
%  ------------------------------------------------------------------

   if (INTCPT == 1)
     X = [ ones(n,1) X ]   ; 
   end
 
    XTXI = inv(X' * X);
    Coefficients = XTXI * X' * Y ;

%  Calculate the fitted regression values
%  --------------------------------------

    Y_hat = X * Coefficients;

%  Calculate R-squared
%  -------------------
% The calculation used for R-squared and the F-statistic here are based
% on the total, un-corrected sum of the squares as describe by Neter and
% Myers. Note that the meaning of R-squared changes for the case of
% regression without a y-intercept. This approach will yield the same
% results as SysStat, SAS, SPSS and BMDP but will differ from that of
% Excel, Quattro Pro, and the MATLAB regress.m function (for the case of
% no y-intercept in the model -- all packages including this one will
% agree for the case of linear regression with a
% y-intercept). Essentially, it is wise to find a way to
% keep the y-intercept (even if it is near zero) in the model to analyze
% it in a meaningful way that everyone can understand.

if (INTCPT == 1)

   RSS = norm(Y_hat - mean(Y))^2;   % Regression sum of squares.
   TSS = norm(Y - mean(Y))^2;       % Total sum of squares (regression plus residual).
   R_sq = RSS / TSS;                % R-squared statistic.

else

   RSS = norm(Y_hat)^2;             % Regression sum of squares.
   TSS = norm(Y)^2;                 % Total, un-corrected sum of squares.
   R_sq = RSS / TSS;                % R-squared statistic.

end

% $$$ % Alternative calculation of R-squared
% $$$ % ====================================
% $$$ % The follwing equation is from Judge G, et al. "An Introduction to the theory
% $$$ % and practice of econometrics", New York : Wiley, 1982. It is the
% $$$ % squared (Pearson) correlation coefficient between the predicted and
% $$$ % dependent variables. It is the same equation regardless of whether an
% $$$ % intercept is included in the model; however, it may yield a negative
% $$$ % R-squared for a particularily bad fit.
% $$$ covariance_Y_hat_and_Y = (Y_hat - mean(Y_hat))' * (Y - mean(Y));
% $$$ covariance_Y_hat_and_Y_hat = (Y_hat - mean(Y_hat))' * (Y_hat - mean(Y_hat));
% $$$ covariance_Y_and_Y = (Y - mean(Y))' * (Y - mean(Y));
% $$$ R_sq = (covariance_Y_hat_and_Y / covariance_Y_hat_and_Y_hat) * ...
% $$$        (covariance_Y_hat_and_Y / covariance_Y_and_Y);


%  Calculate residuals and standard error
%  --------------------------------------

    residuals = Y - Y_hat;

    if (INTCPT == 1)
       S_err = sqrt(residuals' * residuals / (n - k - 1) );
    else
       S_err = sqrt(residuals' * residuals / (n - k) );
    end

%  Calculate the standard deviation and t-values for the regression coefficients
%  -----------------------------------------------------------------------------

    covariance = XTXI .* S_err^2;
    
    C = sqrt(diag(covariance, 0));

    % (n.b. Need to perform a 2-tailed t-test)
    % ****************************************
    if (INTCPT == 1)
      p_value = 2 * (1 - tcdf(abs(Coefficients./C), (n - (k + 1))));
    else
      p_value = 2 * (1 - tcdf(abs(Coefficients./C), (n - k)));
    end

    Coef_stats = [ Coefficients, C, (Coefficients./C), p_value];


% Estimator of error variance.
% ----------------------------

if (INTCPT == 1) 

     SSR_residuals = norm(Y - Y_hat)^2;
     TSS = norm(Y - mean(Y))^2;     % Total sum of squares (regression plus residual).

     F_val = (TSS - SSR_residuals) / k / ( SSR_residuals / (n - (k + 1)));

     F_val = [F_val (1 - fcdf(F_val, k, (n - (k + 1)))) ];

else

     SSR_residuals = norm(Y - Y_hat)^2;
     TSS = norm(Y)^2;                % Total sum of squares (regression plus residual).

     F_val = (TSS - SSR_residuals) / k / ( SSR_residuals / (n-k));

     F_val = [F_val (1 - fcdf(F_val, k, (n - k))) ];

end


