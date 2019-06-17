

## Common functions 

```@autodocs
Modules = [DigitalComm]
Pages   = ["DigitalComm.jl"]
Order   = [:function, :type]
Depth	= 1
```


## Quadrature Amplitude Modulation 

```@autodocs
Modules = [DigitalComm]
Pages   = ["genBitSequence.jl","bitMapping.jl","bitDeMapping.jl","hardConstellation.jl","symbolDemapper.jl"]
Order   = [:function, :type]
Depth	= 1
```


## Channels 

```@autodocs
Modules = [DigitalComm]
Pages   = ["addNoise.jl"]
Order   = [:function, :type]
Depth	= 1
```

# Waveforms 

## Common functions 

```@autodocs
Modules = [DigitalComm]
Pages   = ["genSig.jl","genZCSequence.jl","getLTEAlloc.jl"]
Order   = [:function, :type]
Depth	= 0
```

## BF-OFDM 

```@autodocs
Modules = [DigitalComm]
Pages   = ["Waveforms/BFOFDM/BFOFDM_filter.jl","Waveforms/BFOFDM/bfofdmSigGen.jl","Waveforms/BFOFDM/bfofdmSigDecode.jl","Waveforms/BFOFDM/carrierManipulation.jl"]
Order   = [:function, :type]
Depth	= 0
```

## FBMC 

```@autodocs
Modules = [DigitalComm]
Pages   = ["Waveforms/FBMC/fbmcSigGen.jl","Waveforms/FBMC/fbmcSigDecode.jl"]
Order   = [:function, :type]
Depth	= 0
```


## OFDM 

```@autodocs
Modules = [DigitalComm]
Pages   = ["Waveforms/OFDM/ofdmSigGen.jl","Waveforms/OFDM/ofdmSigDecode.jl"]
Order   = [:function, :type]
Depth	= 0
```

## SC-FDMA 

```@autodocs
Modules = [DigitalComm]
Pages   = ["Waveforms/UFOFDM/ufofdmSigGen.jl","Waveforms/UFOFDM/ufofdmSigDecode.jl"]
Order   = [:function, :type]
Depth	= 0
```

##  UF-OFDM

```@autodocs
Modules = [DigitalComm]
Pages   = ["Waveforms/filterUFOFDM.jl","Waveforms/SCFDMA/scfdmaSigGen.jl","Waveforms/SCFDMA/scfdmaSigDecode.jl"]
Order   = [:function, :type]
Depth	= 0
```


## WOLA 

```@autodocs
Modules = [DigitalComm]
Pages   = ["Waveforms/WOLA/wolaSigGen.jl","Waveforms/WOLA/wolaSigDecode.jl"]
Order   = [:function, :type]
Depth	= 0
```
