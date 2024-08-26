# SEAMOUR Simulink Model Simulation

This repository contains the MATLAB code to simulate the SEAMOUR robot, a bio-robotic sea lion, in a Simscape environment. The simulation includes four different strokes acting on the system:

1. **Characteristic Sea Lion Stroke**: This stroke is based on the natural swimming motion of an actual sea lion.
2. **Learned Stroke A**: Developed using a Soft Actor-Critic (SAC) algorithm.
3. **Learned Stroke B**: Another variant developed using the SAC algorithm.
4. **Learned Stroke C**: A third variant developed using the SAC algorithm.

The results of these simulations are presented in the paper titled **"Using Reinforcement Learning to Develop a Novel Gait for a Bio-robotic California Sea Lion"**, co-authored by Anthony C. Drago, Shraman Kadapa, Nicholas Marcouiller, Harry Kwatny, and James Tangorra. The paper was published in the *Biomimetics* journal.

## How to Run the Simulation

1. Ensure you have MATLAB and Simscape installed.
2. Clone this repository to your local machine.
3. Open MATLAB and navigate to the directory where you cloned the repository.
4. Run the `SEAMOUR_Sim` Simulink model by executing the script 'Test_Sealion_Gaits'. This will simulate the different strokes on the SEAMOUR robot.
5. This is just meant to display the outcomes of these 4 strokes.  The learning environment is separate.

## Dependencies

To run this code, you will need the following MATLAB toolboxes and packages:

- Matlab 2023a
- Simulink
- Simscape
- Optimization Toolbox
- Signal Processing Toolbox
- Curve Fitting Toolbox

## Usage

After running the simulation, the results will be plotted and displayed, showing the effects of each stroke on the SEAMOUR robot's motion. You can compare the performance of the natural sea lion stroke against the three strokes developed using reinforcement learning.
