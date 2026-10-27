# Shiny Apps & Interactive Tools

This repository contains a collection of interactive tools and applications developed using **R and Shiny**. The tools are designed to make statistical analyses, calculations, and other workflows more accessible through user-friendly web interfaces.

## Applications

The repository currently contains four integrated tools:

### 1. BayesApp
An interactive application that illustrates that sensitivity and specificity might vary across pre-test probabilities.

### 2. UtilityApp
An interactive application that estimates the clinical usefulness of additional testing, based on AUC or sensitivity/specificity in a user-defined context.

### 3. MP_iNBApp
An interactive application that estimates the potential value in terms of individualized Net Benefit of additional MammaPrint testing in addition to Breast PREDICT. It aims to assist in the discussion of testing decisions in general.

### 4. iNBApp
An interactive application that estimates the Net Benefit and individualized Net Benefit of two competing test in a user-defined context.

## Accessing the Apps

The applications are developed in RStudio, version-controlled through GitHub, and published online using Posit Connect Cloud.

- **BayesApp:** [Open application](https://01a17df9-ae31-12df-1f65-29a8ed934cbf.share.connect.posit.cloud/)
- **UtilityApp:** [Open application](https://01a17e28-5f6b-bfdd-9335-71df7754417d.share.connect.posit.cloud/)
- **MP_iNBApp:** [Open application](https://01a17e35-a930-9905-4dc2-dcb2ddde3d1f.share.connect.posit.cloud/)
- **iNBApp:** [Open application](https://01a17e2f-773e-3064-08d0-85610554076a.share.connect.posit.cloud/)

## Repository Structure

Each application is maintained in its own subdirectory.

```text
Code-Shinyapps/
├── BayesApp/
├── UtilityApp/
├── MP_iNBApp/
├── iNBApp/
├── Code-Shinyapps.Rproj
└── README.md
```

Each app folder contains the files needed to run its respective application.

## Running an Application Locally

1. Download the relevant files.
2. Open the relevant application folder in RStudio.
3. Install any required R packages that are not already installed.
4. Run the app using RStudio's **Run App** button or execute:

   ```r
   shiny::runApp("app.R")
   ```

   Run this command from the relevant application directory.

## Technologies

- [R](https://www.r-project.org/)
- [Shiny](https://shiny.posit.co/)
- [RStudio](https://posit.co/products/open-source/rstudio/)
- [GitHub](https://github.com/)
- [Posit Connect Cloud](https://connect.posit.cloud/)

## Development and Updates

Applications are maintained in this repository using Git and GitHub. When updates are done they are subsequently updated to Posit Connect Cloud.

## Disclaimer

The applications in this repository are provided for informational, educational, and exploratory purposes only. It must not be used to support decision making, provide conclusive evidence of usefulness or to provide medical or diagnostic services. Although care has been taken in their development, no guarantee is made regarding the accuracy, completeness, reliability, or suitability of the results.

