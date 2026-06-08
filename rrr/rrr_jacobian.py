import numpy as np

A1 = 0.4
A2 = 0.3
A3 = 0.2

DT = 0.05

def sysCall_init():
    sim = require('sim')

    self.base_frame = sim.getObject("/base")

    self.joint_hdl = []
    self.joint_hdl.append(sim.getObject("/j1"))
    self.joint_hdl.append(sim.getObject("/j2"))
    self.joint_hdl.append(sim.getObject("/j3"))

    self.X_d_hdl = sim.getObject("/goal")
    self.X_d = sim.getObjectPose(self.X_d_hdl, self.base_frame)

    self.joint = []
    self.joint.append(sim.getJointPosition(self.joint_hdl[0]))
    self.joint.append(sim.getJointPosition(self.joint_hdl[1]))
    self.joint.append(sim.getJointPosition(self.joint_hdl[2]))

def sysCall_sensing():
    self.joint[0] = sim.getJointPosition(self.joint_hdl[0])
    self.joint[1] = sim.getJointPosition(self.joint_hdl[1])
    self.joint[2] = sim.getJointPosition(self.joint_hdl[2])

    self.X_d = sim.getObjectPose(self.X_d_hdl)

def sysCall_actuation():
    # pegamos o valor de onde esta o efetuador agora, o equivalente a fazer cinematica direta
    X_c = sim.getObjectPose(sim.getObject("/endpoint"))

    # tratamento para o espaco de trabalho util do robo (com conversao para quaternio)
    x_c = np.array([X_c[0], X_c[1], 2*np.arctan2(X_c[5], X_c[6])])
    x_d = np.array([self.X_d[0], self.X_d[1], 2*np.arctan2(self.X_d[5], self.X_d[6])])

    # calculamos a diferenca
    Dx = x_d - x_c
    print(f"Dx : {Dx}")

    dx = Dx if np.linalg.norm(Dx) < 0.1 else (Dx/np.linalg.norm(Dx)) * 0.1
    print(f"dx : {dx}")

    dx = dx.reshape((3,1))

    # Calcular inversa da jacobiana
    J = jacobian(self.joint)
    print(f"J: {J}")
    print(f"det(J): {np.linalg.det(J)}")

    #J_inv = np.linalg.inv(J)
    #J_inv = np.linalg.pinv(J)
    J_inv = np.linalg.inv(J.T @ J + 0.01*np.eye(3))@J.T

    dq = J_inv @ dx

    dq = dq.ravel()

    q = np.array(self.joint) + dq*DT

    # Setar o valor das juntas do robo
    for i, theta in enumerate(q):
        sim.setJointTargetPosition(self.joint_hdl[i], theta)


def sysCall_cleanup():
    # do some clean-up here
    pass

def jacobian(q):
    J = np.array([[-A1*np.sin(q[0]) -A2*np.sin(q[0]+q[1]) -A3*np.sin(q[0]+q[1]+q[2]), -A2*np.sin(q[0]+q[1]) -A3*np.sin(q[0]+q[1]+q[2]), -A3*np.sin(q[0]+q[1]+q[2])],
        [A1*np.cos(q[0]) +A2*np.cos(q[0]+q[1]) +A3*np.cos(q[0]+q[1]+q[2]), +A2*np.cos(q[0]+q[1]) +A3*np.cos(q[0]+q[1]+q[2]), A3*np.cos(q[0]+q[1]+q[2])],
        [1, 1, 1]])
    return J

# See the user manual or the available code snippets for additional callback functions and details
