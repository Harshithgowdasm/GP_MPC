# Data Driven Model Predictive Control using Gaussian Processes

### Project Description:
 Model predictive control (MPC) has gained more attention in the past years as it
 is a powerful modern control technique that relies on repeatedly solving an open-loop
 optimal control problem. In order to obtain high-quality optimal input signals, an accurate
 prediction model is required, which is often unavailable or difficult to identify. In recent
 years data-driven control approaches have become more interesting because they use
 only measured data to control unknown systems without prior model identification. A
 major issue in data-driven control is the treatment of measurement noise, which in recent
 publications is treated to be deterministically bounded. However, this assumption
 does not generally hold for a real system. Alternatively, a stochastic model can be used
 to model the state transition function. However, there is a lack of theoretical
 guarantees of probabilistic models, e.g., Gaussian processes.

 The goal of this thesis is to implement a probabilistic MPC scheme relying on Gaussian
 process models to represent the transition function of a discrete-time system with noisy
 measurements and uncertainty models. The thesis proves that one probabilistic GP-MPC
 algorithm is enough for both linear and nonlinear models. A common challenge is that most
 data-driven approaches require numerous interactions with the environment. However,
 conducting a vast number of interactions can often be unfeasible in real-world scenarios,
 especially in fields like robotics. The proposed GP-MPC model is capable of controlling
 any system with a small number of interactions, which increases the data efficiency even
 in an uncertainty model.


### Approach
 A promising approach to enhance the data efficiency of the Model-Based Learning Model
 without relying on task-specific prior knowledge is to learn models of the underlying system
 dynamics. When a high-quality model is available, it can effectively serve as a stand-in for
 the real environment. This means that beneficial policies can be derived from the model
 without the need for additional interactions with the actual system. However, accurately
 modeling the underlying transition dynamics presents a significant challenge and inevitably
 introduces model errors.

 To address these model errors, probabilistic models have been proposed specifically the
 Probabilistic Gaussian process. These models explicitly incorporate uncertainties about
 the system dynamics. By considering model uncertainty, Model-based learning algorithms
 can substantially reduce the number of interactions required with the real system. This
 approach acknowledges that while perfect models are often unattainable, probabilistic
 models offer a more robust framework for navigating the inherent uncertainties in the
 modeling process. Thus, leveraging probabilistic models can lead to more efficient learning
 and improved performance in learning tasks by better accommodating the inevitable
 inconsistencies between the learned model and the real system.

 This probabilistic Gaussian process has a limitation in that it cannot accept Gaussian input
 and predict over the entire prediction horizon (detailed explanation in Section 2.3). To
 address this, a transition model is proposed for long-term prediction, which serves as input
 to the MPC controller. However, an open-loop controller cannot stabilize the system on
 its own. Therefore, obtaining a feedback controller is essential. Model Predictive Control
 (MPC) provides a practical framework for this purpose. During interaction with the
 system, MPC determines an N-step open-loop control trajectory u0,...,uN−1, starting
 from the current state xt. Only the first control signal u0 is applied to the system.
 When the system transitions to xt+1, the Gaussian Process (GP) model is updated with
 the newly available information, and MPC re-plans u0,...,uN−1. This process effectively
 transforms an open-loop controller into an implicit closed-loop (feedback) controller by
 continuously re-planning N steps ahead from the current state, as illustrated in Figure: 
    <p align="middle">
        <img src="https://github.com/Harshithgowdasm/PA_Harshith_Gowda/blob/main/papers/blockdia.png" width="80%" />
    </p>


### Examples
The application of GP-MPC framework lies in its versatility, as it can be applied to a wide range of systems,
 whether linear or nonlinear, with or without uncertainty, and with or without noisy
 measurements. Model uncertainty refers to parameter variations and external disturbances
 that affect the system’s behavior. These mathematical uncertainty models represent all
 aspects of real-world systems, leading to uncertainty in predicting their behavior. To
 demonstrate its effectiveness, we have chosen to implement the GP-MPC framework on
 both a linear system, the DC motor, and a nonlinear system, the Van der Pol oscillator.

  ### Linear system - DC Motor

      <p align="middle">
        <img src="https://github.com/Harshithgowdasm/PA_Harshith_Gowda/blob/main/examples/Dc%20motor%20plant/Results/Gp_MPC_without_uncertanity/gp_mpc_wo_1.svg" width="80%" />
    </p>
  