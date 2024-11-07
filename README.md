[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Model Averaging under Flexible Loss Functions

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE.txt).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
[Model Averaging under Flexible Loss Functions](https://doi.org/10.1287/ijoc.2023.0291) by D. Gu, Q. Liu, and X. Zhang. 

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0291

https://doi.org/10.1287/ijoc.2023.0291.cd

Below is the BibTex for citing this snapshot of the repository.

```
@article{gu2024model,
  author =        {D. Gu, Q. Liu, and X. Zhang},
  publisher =     {INFORMS Journal on Computing},
  title =         {{Model Averaging under Flexible Loss Functions}},
  year =          {2024},
  doi =           {10.1287/ijoc.2023.0291.cd},
  url =           {https://github.com/INFORMSJoC/2023.0291},
  note =          {Available for download at https://github.com/INFORMSJoC/2023.0291},
}  
```

## Description

To address model uncertainty under flexible loss functions in prediction problems, we propose a model averaging method that accommodates various loss functions,
 including asymmetric linear and quadratic loss functions, as well as many other asymmetric/symmetric loss functions as special cases. The flexible loss function allows the proposed method to average a large range of models, such as the quantile and expectile regression models. To determine the weights of the candidate models, we establish a J-fold cross-validation criterion. Asymptotic optimality and weights convergence are proved for the proposed method. Simulations and an empirical application show the superior performance of the proposed method, compared with other methods of model selection and averaging.

This project contains three folders: `data`, `results`, `src`.

- `data`: This folder includes the data used in the paper.
- `src`: This folder contains the source code for the simulations and two empirical applications.
- `results`: This folder contains the results of the experiments.

## Results

All detailed results are available in the `results` folder.

## Data

A description of the datasets used for the empirical study can be found in the `data` folder.

## Replicating

To reproduce each result in the paper, please run corresponding file in `src` . For example, to perform the simulations of DGP 1 in Simulation Design I, please run R

`src/DGP1.R`

See the README.md file in `src` for a detailed description.

Please note that we utilized large-scale computing equipment for multi-core parallel processing. For standard computers without these capabilities, `sfLapply` can be replaced by `lapply` for smaller-scale calculations.

