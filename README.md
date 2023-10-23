# adaptive-NN-control-manipulator-in-task-space
Referring to the control methods in [1] and [2], we study the design method of adaptive control based on RBF (Radial Basis Function) network for dual-joint manipulator in task space. The approach eliminates the need for estimation of the dynamic model and avoids time-consuming training process. By introducing the GL (Ge-Lee) matrix and its multiplication operator, the control law is obtained using direct parameter identification without the need for the inverse of Jacobian matrix. The approach employs robust control to suppress modeling error and bound perturbation.
reference: 
[1] Liu JinKun. Robot Control System Design and MATLAB Simulation[M]. Tsinghua University Press, 2008.
[2] Shuzhi S G, Hang C C, Woon L C. Adaptive neural network control of robot manipulators in task space[J]. IEEE transactions on industrial electronics, 1997, 44(6): 746-752.