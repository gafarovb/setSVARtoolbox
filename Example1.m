
% Exmaple 1;
% standard dataset  - MSG

msgModel = SVAR;
estimatedBounds = msgModel.onesidedLowerIRFHatAnalytic;

ConfidenceLevel = 0.95;
tic;
CSanalytic = msgModel.onesidedUpperIRFCSAnalytic( ConfidenceLevel);
toc;

tic;
CSbonferroni = msgModel.twosidedIRFCSbonferroni( ConfidenceLevel);
toc;