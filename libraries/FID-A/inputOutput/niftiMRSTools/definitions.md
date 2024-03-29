# Documentation for the definitions.json file

The `definitions.json` file contains machine-readable definitions of NIfTI-MRS standard meta data and dimension tags.

Type fields should be generic JSON types: `number`, `string`, `array`, `object`, `bool`
or an array indicating a (optionally nested) array type and element type : `[array, number]` for a vector of numbers or `[array, array, string]` for a matrix of strings.

## nifti-mrs version number.
Two numeric elements, major and minor

## Possible dimension tags and descriptions
An object with keys as the possible tage, with documentation/use of those tags as the value.

## Header fields
All header fields are defined with a JSON object containing type information `type`, appropriate unit definitions `units`, a dco string `doc`, and whether the field should be removed during anonymisation `anon`.

### Required
Two fields are required, `SpectrometerFrequency` and `ResonantNucleus`.

### Standard defined
These fields are optional but must not be redefined.
This section currently contains:
#### 5.1 MRS specific Tags
- SpectralWidth
- EchoTime
- RepetitionTime
- InversionTime
- MixingTime
- AcquisitionStartTime
- ExcitationFlipAngle
- TxOffset
- VOI
- WaterSuppressed
- WaterSuppressionType
- SequenceTriggered
#### 5.2 Scanner information
- Manufacturer
- ManufacturersModelName
- DeviceSerialNumber
- SoftwareVersions
- InstitutionName
- InstitutionAddress
- TxCoil
- RxCoil
#### 5.3 Sequence information
- SequenceName
- ProtocolName
#### 5.4 Sequence information
- PatientPosition
- PatientName
- PatientID
- PatientWeight
- PatientDoB
- PatientSex
#### 5.5 Provenance and conversion metadata
- ConversionMethod
- ConversionTime
- OriginalFile
#### 5.6 Spatial information
- kSpace
#### 5.7 Editing Pulse information structure
- EditCondition
- EditPulse
#### 5.8 Processing Provenance
- ProcessingApplied
