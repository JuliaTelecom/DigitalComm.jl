""" 
---  
Apply Chebyshev polynom of order L to input x
# --- Syntax 
	y = cheb(L,x)
# --- Input parameters 
- L	  : Chebyshev order [Int] 
- x	  : Input [Any]
# --- Output parameters 
- y	  : Output 
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function cheb(L,x)
	## Apply Chebyshev polynom of order m to input x
	 x[abs.(x).<=1] = cos.(L .* acos.(x[abs.(x).<=1]));
	 x[abs.(x).>1]  = cosh.(L .* acosh.(abs.(x[abs.(x).>1])));
	 return x;
end

""" 
---  
Returns the Dolp Chebyshev filter of order L with desired attenation att

See Peter Lynch, "The Dolph-Chebyshev Window: A Simple Optimal Filter", Monthly Weather Review, Vol. 125, pp. 655-660, April 1997. (http://www.maths.tcd.ie/~plynch/Publications/Dolph.pdf)
Dolph, "A current distribution for broadside arrays which optimizes the relationship between beam width and side-lobe level", Proc. IEEE, 34, pp. 335-348.
		cheb(m-1, beta * cos(pi * k/m))
W(k) =	-------------------------------
			  cheb(m-1, beta)
# --- Syntax 
 filterDC = dolphChebyshev(n,at)
# --- Input parameters 
- n		  : Size of desired filter 
- at	  : Attenation in dBB 
# --- Output parameters 
- filterDC : Filter impulse response [Array{Float64}] 
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function dolphChebyshev(n,at);
	# --- Filter in frequency domain
	# --- DC constant
    gamma = 10^(-at/20);
    beta  = cosh(1/(n-1) * acosh(1/gamma));
	k     = collect(0:n-1);
    x     = beta*cos.(pi*k/n);
	p	  = cheb(n-1, x);
	if mod(n,2) == 1
		w = real(fft(p));
		M = Int((n+1)/2);
		w = w[1:M]/w[1];
		w = [w[M:-1:2];w];
	else
		off = exp.(1im*pi/n.*collect(0:n-1));
		# half-sample delay (even order)
		p .= p * transpose(off);
		w = real(fft(p));
		M = Int(n/2+1);
		w = w./w[2];
		w .= [w[M:-1:2];w[2:M]];
	end
	return w
end
