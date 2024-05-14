# Service Provider Planning Tool

## Introduction
This repository contains a MATLAB project aimed at simplifying the planning process for a service provider owning 340 channels in the 900 MHz band.

## Requirements
- MATLAB installed on your system.
- Basic understanding of cellular network planning concepts.

## Usage

- **Clone the Repository:**
    ```bash
    git clone [https://github.com/yourusername/planning-tool.git](https://github.com/Ahmed-Salah-312/Wireless-Communication-Networks)
    ```

- **Input Parameters:**
    - Follow the on-screen prompts to input parameters like Grade of Service (GOS), city area, user density, SIRmin, and sectorization method.

- **View Results:**
    - Check the generated design parameters and plots.

## Design Parameters

- **Cluster Size**: 
- **Number of Cells**: 
- **Cell Radius**: 
- **Traffic Intensity per Cell and per Sector**: 
- **Base Station Transmitted Power**: 
- **Received Power Plot**: 

## Validation Figures

### Figure 1: Cluster Size vs. SIRmin
- Curves for omni-directional, 120° sectorization, and 60° sectorization designs.

### Figure 2: 
- Plot the number of cells versus GOS (1% to 30%) for SIRmin = 19dB & user density = 1400 users/km^2.
- Plot the traffic intensity per cell versus GOS (1% to 30%) for SIRmin = 19dB & user density = 1400 users/km^2.

### Figure 3: 
- Plot the number of cells versus GOS (1% to 30%) for SIRmin = 14dB & user density = 1400 users/km^2.
- Plot the traffic intensity per cell versus GOS (1% to 30%) for SIRmin = 14dB & user density = 1400 users/km^2.

### Figure 4: 
- Plot the number of cells versus user density (100 to 2000 users/km^2) for SIRmin = 14dB & GOS = 2%.
- Plot the cell radius versus user density (100 to 2000 users/km^2) for SIRmin = 14dB & GOS = 2%.

### Figure 5: 
- Plot the number of cells versus user density (100 to 2000 users/km^2) for SIRmin = 19dB & GOS = 2%.
- Plot the cell radius versus user density (100 to 2000 users/km^2) for SIRmin = 19dB & GOS = 2%.


## Acknowledgments

- Thanks to the contributors who helped in developing this planning tool.

## Disclaimer

This tool is provided as-is and without any warranties. Use at your own risk.
