# mec-with-ris-control

Code implementing the performance evaluation of the paper "Control Aspects for Using RIS in Latency-Constrained Mobile Edge Computing"

## Abstract
> This paper investigates the role and the impact of control operations for dynamic mobile edge computing (MEC) empowered by Reconfigurable Intelligent Surfaces (RISs), in which multiple devices offload their computation tasks to an access point (AP) equipped with an edge server (ES), with the help of the RIS. While usually ignored, the control aspects related to channel estimation (CE), resource allocation (RA), and control signaling play a fundamental role in the user-perceived delay and energy consumption. In general, the higher the resources involved in the control operations, the higher their reliability; however, this introduces an overhead, which reduces the number of resources available for computation offloading, possibly increasing the overall latency experienced. Conversely, a lower control overhead translates to more resources available for computation offloading but impacts the CE accuracy and RA flexibility. This paper establishes a basic framework for integrating the impact of control operations in the performance evaluation of the RIS-aided MEC paradigm, clarifying their trade-offs through theoretical analysis and numerical simulations.

The paper has been presented to Asilomar Conference on Signals, Systems, and Computers 2023, Asilomar, Pacific Grove, CA, USA.


Channel generation needs to be performed by running
```
channel_generation.m
```

The main performance are obtainable by running
```
main_ce.m
main_error_proba.m
main_slot_time.m
```
Note that the scripts save `.mat` files in specific directories. Create the directory **before running the scripts**.

