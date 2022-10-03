# windyPlume 

windyPlume is a collection of Python scripts that can be used to calibrate images
of wind-affected volcanic plumes. For a full description of how the software
works see a) comments in code and b) manuscript by Snee et al.

## Contents:

Script_2_Calibrate.py - This is the primary script, to be ran by the user.
                        Variables in the User-Input section are to be edited by
			the user.

calibrate_function.py - Function that performs the geometric calibration based on
                        Bombrun et al. (2018)

Calibrate_allTogether_function.py - Function that performs the calibration

calcInclination_function.py - Function to estimate the camera inclination (if it
                              is not known)

calibrateWind_function.py - Function to perform the wind correction to the
                            calibration

extractWeather_function.py - Function to extract the weather data from netCDF
                             files

getAllUncertainty_function.py - Function called by Script_2_Calibrate.py -
                                Performs the calubration for different values of
				input parameters corresponding to the ranges of
				their uncertainties

## Pre-requisites

windyPlume requires Python3 as well as the following modules:

1. numpy
2. matplotlib
3. scipy
4. math
5. netCDF4
6. datetime
7. sys
8. copy

These modules can be installed using pip:

$ pip install <packageName>

e.g.,

$ pip install numpy

## Issues

Reporting of issues and/or contacting the authors can be achieved through the
GitHub repository for the code:
https://github.com/Paul-Jarvis/windyPlume/issues - THE CORRECT URL NEEDS TO BE ENTERED HERE. To report an issue, click on
the "New issue" icon and leaving a comment. If reporting a bug or unexpected
behaviour, please provide as much detail as possible, inlcuding a minimal
example. The authors will make a best endeavour to respond in a reasonable time

## Contributing to the software

Users are welcome to edit the software for their own purposes and are welcome to
contact the authors during this process (see "Issues" above). We invite such
users to contribute their edits to our code. The easiest way to achieve this is
probably to:

1. Create a fork of the repository. To do this, go to the code repository at
https://github.com/Paul-Jarvis/windyPlume (NOT CORRECT URL)and click on the "Fork" icon in the
top right. More information about how to create a fork and work with the
subsequent copy of the repository can be found at
https://docs.github.com/en/get-started/quickstart/fork-a-repo.

2. Work with your fork to make the desired changes.

3. Submit a pull request, proposing your changes, at
https://github.com/Paul-Jarvis/windyPlume/pulls (NOT CORRECT URL) by clicking on the "New pull
request" icon. More information on pull requests can be found at
https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork.

## Authors

Eveanjelene Snee
Paul Jarvis

## References

Snee, et al. (in prep) - Image analysis of volcanic plumes: A simple calibration
tool to correct for the effect of wind

## License - This project is licensed under the MIT license