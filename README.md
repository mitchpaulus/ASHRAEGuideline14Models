# ASHRAEGuideline14Models
This repository contains different linear change points related to ASHRAE Guideline 14. 

## Available Models 

The three parameter cooling and heating models are currently available.

```Matlab
[coefficients, minSSE] = threeparametercooling(x, y);
[coefficients, minSSE] = threeparameterheating(x, y);
```

The four and five parameter models are currently not available.
