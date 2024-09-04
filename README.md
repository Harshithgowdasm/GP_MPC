# Data Driven Model Predictive Control using Gaussian Processes

### Project Description:
 Model predictive control (MPC) has gained more attention in the past years as it is a powerful modern control technique that relies on repeatedly solving an open-loop optimal control problem. In order to obtain high-quality optimal input signals, an accurate prediction model is required, which is often unavailable or difficult to identify. In recent years data-driven control approaches have become more interesting because they use only measured data to control unknown systems without prior model identification. A major issue in data-driven control is the treatment of measurement noise, which in recent publications is treated to be deterministically bounded. However, this assumption does not generally hold for a real system. Alternatively, a stochastic model can be used to model the state transition function. However, there is a lack of theoretical guarantees of probabilistic models, e.g., Gaussian processes.

 The goal of this thesis is to implement a probabilistic MPC scheme relying on Gaussian process models to represent the transition function of a discrete-time system with noisy measurements and uncertainty model The thesis proves that one probabilistic GP-MPC algorithm is enough for both linear and nonlinear models. A common challenge is that most data-driven approaches require numerous interactions with the environment. However, conducting a vast number of interactions can often be unfeasible in real-world scenarios, especially in fields like robotics. The proposed GP-MPC model is capable of controlling any system with a small number of interactions, which increases the data efficiency evenin an uncertainty model.


### Approach
A promising approach to enhance the data efficiency of Model-Based Learning without relying on task-specific prior knowledge is to learn models of the underlying system dynamics. A high-quality model can act as a stand-in for the real environment, enabling policy derivation with reduced interactions with the actual system. However, accurately modeling transition dynamics introduces challenges and potential errors.

To address these errors, probabilistic models like the Gaussian Process (GP) have been proposed, incorporating uncertainties in system dynamics. This reduces the number of required interactions by accounting for model uncertainty. Although GP models are robust, they struggle with Gaussian inputs and long-term predictions. To overcome this, a transition model is introduced for long-term prediction, feeding into a Model Predictive Control (MPC) framework. MPC generates an N-step open-loop control trajectory, applying only the first control signal and updating the GP model with new information as the system evolves. This process transforms the open-loop controller into an implicit feedback controller through continuous re-planning.This process effectively transforms an open-loop controller into an implicit closed-loop (feedback) controller by continuously re-planning N steps ahead from the current state, as illustrated in Figure: 
    <p align="middle">
        <img src="https://github.com/Harshithgowdasm/PA_Harshith_Gowda/blob/main/papers/blockdia.png" width="80%" />
    </p>


### Examples
The application of GP-MPC framework lies in its versatility, as it can be applied to a wide range of systems, whether linear or nonlinear, with or without uncertainty, and with or without noisy measurements. Model uncertainty refers to parameter variations and external disturbances that affect the system’s behavior. These mathematical uncertainty models represent all aspects of real-world systems, leading to uncertainty in predicting their behavior. To demonstrate its effectiveness, we have chosen to implement the GP-MPC framework on both a linear system, the DC motor, and a nonlinear system, the Van der Pol oscillator.

  ### Linear system - DC Motor
    
 ### GP-MPC without noise   
<p align="middle">
        <img src="https://github.com/Harshithgowdasm/PA_Harshith_Gowda/blob/main/examples/Dc%20motor%20plant/Results/Gp_MPC_without_uncertanity/gp_mpc_wo_1.svg" width="80%" />
</p>

 ### GP-MPC with noise
<p align="middle">
        <img src="https://github.com/Harshithgowdasm/PA_Harshith_Gowda/blob/main/examples/Dc%20motor%20plant/Results/GP_MPC_with_uncertainity/mpc_dc_noise_1.svg" width="80%" />
</p>

### Non Linear system - Van der Pol Oscillator
 ### GP-MPC without noise   
<p align="middle">
        <img src="https://github.com/Harshithgowdasm/PA_Harshith_Gowda/blob/main/examples/Van_der_pol_oscillator/results/GP_MPC_without_uncertanity/GP_mpc_better.svg" width="80%" />
</p>

 ### GP-MPC with noise
<p align="middle">
        <img src="https://github.com/Harshithgowdasm/PA_Harshith_Gowda/blob/main/examples/Van_der_pol_oscillator/results/GP_MPC_with_uncertanity/2_Gp_mpc_uncert_vdp.svg" width="80%" />
</p>

 ### Reference tracking 
<p align="middle">
        <img src="https://github.com/Harshithgowdasm/PA_Harshith_Gowda/blob/main/examples/Van_der_pol_oscillator/results/GP_MPC_with_uncertanity/4_gp_mpc_set_1_0.svg " width="80%" />
</p>

 