# Polyomino self-assembly core repository

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/8eae273ab3f34d77b3be158a33613ce6)](https://app.codacy.com/app/ASLeonard/polyomino_core?utm_source=github.com&utm_medium=referral&utm_content=ASLeonard/polyomino_core&utm_campaign=Badge_Grade_Dashboard)

Base for polyomino model development

## Getting Started
git clone --recurse-submodules

Some minor diffulties on submodule cloning, but include the flag '--recurse-submodules' when cloning an implementing repository. 

Alternatively, just clone this repository automatically inside an implementing repository

### Extending the polyomino core model
Including this repository as a submodule allows for rapid prototyping of the polyomino model.

Generic assembly can be quickly implemented using existing functions. Extend the model as
  - extend base assembly 
```cpp
class NewAssemblyModel : public PolyominoAssembly<NewAssemblyModel>
```
  - determine the interface type (e.g. integer, string, something more complex, etc.)
```cpp
using Genotype = std::vector<subunit type>
```
  - implement the interaction matrix that returns the interaction strength between two interfaces
```cpp
double NewAssemblyModel::InteractionMatrix(subunit A, subunit B)
```
  - implement any specific methods
    - mutations
    - random initilising
    
### Using the core model
The main calls to the core model are to the assembly and classification
  - get interaction edges
  ```cpp
  const auto edges = NewAssemblyModel::GetActiveInterfaces(genotype);
  ```
  - assemble polyomino (and can track edges)
  ```cpp
  auto assembly_information=NewAssemblyModel::AssemblePolyomino(edges,interacting_indices);
  ```
  - get phenotype from the raw assembly
  ```cpp
  Phenotype phen=GetPhenotypeFromGrid(assembly_information);
  ```
  
  - potentially classify phenotype using a runtime filled table
  
  ```cpp
  Phenotype_ID pid = phenotypeTable->GetPhenotypeID(phen);
  ```
  
More robust examples available in the implementation of `polyomino_interface' model

## Building
The core library is header only, so need to add an include path to whatever building process used
  - -Ipath/to/library/polyomino_core
