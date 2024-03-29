# Function overview

## Common functions

```@autodocs
Modules = [DigitalComm]
Pages   = ["DigitalComm.jl"]
Order   = [:function, :type]
```

## NRZI Encoding

```@autodocs
Modules = [DigitalComm]
Pages   = ["NRZI.jl"]
Order   = [:function, :type]
```

## Quadrature Amplitude Modulation

```@autodocs
Modules = [DigitalComm]
Pages   = ["genBitSequence.jl","bitMapping.jl","bitDeMapping.jl","hardConstellation.jl","symbolDemapper.jl"]
Order   = [:function, :type]
```

## Channels

```@autodocs
Modules = [DigitalComm]
Pages   = ["Channel/addNoise.jl","Channel/rayleighChan.jl","Channel/getChannel.jl", "Channel/addCFO.jl"]
Order   = [:function, :type]
```

## Windows and filters

```@autodocs
Modules = [DigitalComm]
Pages   = ["raisedCosine.jl", "UFOFDM/filterUFOFDM.jl", "WOLA/getWolaWindow.jl"]
Order   = [:function, :type]
```

# Waveforms

## Common functions

```@autodocs
Modules = [DigitalComm]
Pages   = ["genSig.jl","genZCSequence.jl","getLTEAlloc.jl"]
Order   = [:function, :type]
```

## BF-OFDM

```@autodocs
Modules = [DigitalComm]
Pages   = ["Waveforms/BFOFDM/BFOFDM_filter.jl","Waveforms/BFOFDM/bfofdmSigGen.jl","Waveforms/BFOFDM/bfofdmSigDecode.jl","Waveforms/BFOFDM/carrierManipulation.jl"]
Order   = [:function, :type]
```

## FBMC

```@autodocs
Modules = [DigitalComm]
Pages   = ["Waveforms/FBMC/fbmcSigGen.jl","Waveforms/FBMC/fbmcSigDecode.jl"]
Order   = [:function, :type]
```

## OFDM

```@autodocs
Modules = [DigitalComm]
Pages   = ["Waveforms/OFDM/ofdmSigGen.jl","Waveforms/OFDM/ofdmSigDecode.jl"]
Order   = [:function, :type]
```

## SC-FDMA

```@autodocs
Modules = [DigitalComm]
Pages   = ["Waveforms/UFOFDM/ufofdmSigGen.jl","Waveforms/UFOFDM/ufofdmSigDecode.jl"]
Order   = [:function, :type]
```

## UF-OFDM

```@autodocs
Modules = [DigitalComm]
Pages   = ["Waveforms/filterUFOFDM.jl","Waveforms/SCFDMA/scfdmaSigGen.jl","Waveforms/SCFDMA/scfdmaSigDecode.jl"]
Order   = [:function, :type]
```

## WOLA

```@autodocs
Modules = [DigitalComm]
Pages   = ["Waveforms/WOLA/wolaSigGen.jl","Waveforms/WOLA/wolaSigDecode.jl"]
Order   = [:function, :type]
```
