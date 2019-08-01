# Auto_location
This repository contains a series of functions that read in nordic format files that contain P- and S-wave phase picks and perform an automatic location. The automatic location is performed in an iterative manner where in each iteration the picks with the highest residuals are removed until some quality criteria are met (e.g. number of P and S picks, minimum earthquake location rms). To populate the sfiles with phase picks we used the automatic picker called kpick. For the earthquake location we use HYPOCENTER incorporated in SEISAN.

### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```
## Authors

* **Konstantinos Michailos** 
