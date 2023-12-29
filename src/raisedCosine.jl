"""   
Returns the Finite Impulse Response of a Raised Cosine (RC) filter.	The filter is defined by its span (evaluated in number of symbol N), its Roll-Off factor and its oversampling factor. The span corresponds to the number of symbol affected by filter before and after the center point.\n
Output is a Vector{Float64} array of size L= 2KN+1 \n
SRRC definition is based on [1]	 \n
[1]	 3GPP TS 25.104 V6.8.0 (2004-12). http://www.3gpp.org/ftp/Specs/archive/25_series/25.104/25104-680.zip \n
 Syntax \n
		h	= raisedCosine(N,beta,ovS)
Input parameters \n
- N	  	: Symbol span (Int16)
- beta  : Roll-off factor (Float64)
- ovS	: Oversampling rate (Int16)

## Examples
```jldoctest
julia> h = raisedCosine(12,0.5,16);

julia> sum(h)
16.002451134661918
```
"""
function raisedCosine(N,beta,ovS)::Vector{Float64}
	# --- Final size of filter
	nbTaps	= 2 * N * ovS + 1;
	# --- Init output
	h			= zeros(Float64,nbTaps);
	counter		= 0;
	# --- Iterative SRRC definition
	for k = -N*ovS : 1 : N*ovS
		counter		 = counter + 1;
		if k == 0
			## First singular point at t=0
			h[counter]  = 1;#(1-beta) + 3*beta/pi;
		elseif abs(k) == ovS / (2*beta);
			## Second possible singular point
			h[counter]  = pi/4*sin(pi/(2beta))/(pi/(2beta));
		else
			## Classic SRRC formulation (see [1])
			h[counter]	= sin(pi*k/ovS)/(pi*k/ovS) * cos(pi*beta*k/ovS)/(1- 4beta^2*(k/ovS)^2);
		end
	end
	return h
end


"""
Returns the Finite Impulse Response of a Square Root Raised Cosine (SRRC) filter. \n
The filter is defined by its span (evaluated in number of symbol N), its Roll-Off factor and its oversampling factor. The span corresponds to the number of symbol affected by filter before and after the center point.\n
Output is a `Vector{Float64}` array of size L= 2KN+1\n
SRRC definition is based on [1]\n
[1]	 3GPP TS 25.104 V6.8.0 (2004-12). http://www.3gpp.org/ftp/ Specs/archive/25_series/25.104/25104-680.zip\n
Syntax\n
h	= sqrtRaisedCosine(N,beta,ovS) \n
Input parameters \n
- N	  	: Symbol span (Int16)
- beta  : Roll-off factor (Float64)
- ovS	: Oversampling rate (Int16)

## Examples
```jldoctest
julia> h = sqrtRaisedCosine(12,0.5,16);

julia> sum(h)
14.07688841517184
```
"""
function sqrtRaisedCosine(N,beta,ovS)::Vector{Float64}
	# --- Final size of filter
	nbTaps	= 2 * N * ovS + 1;
	# --- Init output
	h			= zeros(Float64,nbTaps);
	counter		= 0;
	# --- Iterative SRRC definition
	for k = -N*ovS : 1 : N*ovS
		counter		 = counter + 1;
		if k == 0
			## First singular point at t=0
			h[counter]  = (1-beta) + 4*beta/pi;
		elseif abs(k) == ovS / (4*beta);
			## Second possible singular point
			h[counter]  = beta/sqrt(2)*( (1+2/pi)sin(pi/(4beta))+(1-2/pi)cos(pi/(4beta)));
		else
			## Classic SRRC formulation (see [1])
			h[counter]  = ( sin(pi*k/ovS*(1-beta)) + 4beta*k/ovS*cos(pi*k/ovS*(1+beta))) / (pi*k/ovS*(1- (4beta*k/ovS)^2) );
		end
	end
    # --- Max @ h[0]
    h   = h ./ maximum(h)
	return h
end
