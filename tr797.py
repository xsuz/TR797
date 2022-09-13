import numpy as np
from numpy import cos,sin,pi,tan
from matplotlib import pyplot as plt


# 定数（既知とするパラメータ）
N = 100 # 分割数
rho = 1.2  # 空気の密度 [kg/m^3]
U = 7.2  # 一様流の流速 [m/s]
Lift=3000.0 # 揚力 [N]
beta=.85
D=10e-2 # 桁の直径
E=3e11 # 桁のヤング率
d=1e-3 # 桁の厚さ
lr=0.1 # 学習率
ITER=500
plt.style.use('dark_background')

class Wing:
    def __init__(self):
        self.l_e=15.0
        self.dS = np.ones(shape=(N))*self.l_e/(N*2) # パネルの半幅
        self.y=np.linspace(self.dS[0],self.l_e-self.dS[N-1],N) # パネルの中心のy座標
        self.z=np.zeros(shape=(N)) # パネルの中心のz座標
        self.phi=np.zeros(shape=(N)) # パネルの上反角
        self.gamma=np.zeros(shape=(N)) # パネル周りの循環
        self.local_bending_moment=np.zeros(shape=(N)) # パネルの中心にかかる曲げモーメント
        self.local_induced_drag=np.zeros(shape=(N)) # パネルに働く誘導抗力
        self.local_lift=np.zeros(shape=(N)) # パネルに働く揚力
        self.V_n=np.zeros(shape=(N)) # パネルの吹きおろし速度
    def optim(self):
        y=self.y
        z=self.z
        dS=self.dS
        l_e=self.l_e
        phi=self.phi

        # Biot-Savartの法則の係数1/2pi*h
        Q=np.zeros(shape=(N,N))
        for i in range(N):
            for j in range(N):
                # パネルj,j'を中心とした回転座標 ... (9)
                yd= (y[i]-y[j])*cos(phi[j])+(z[i]-z[j])*sin(phi[j])
                zd=-(y[i]-y[j])*sin(phi[j])+(z[i]-z[j])*cos(phi[j])
                y2d=(y[i]+y[j])*cos(phi[j])-(z[i]-z[j])*sin(phi[j])
                z2d=(y[i]+y[j])*sin(phi[j])+(z[i]-z[j])*cos(phi[j])
                # パネル間の距離^2 ... (10)
                R_plus2=(yd-dS[j])**2+zd**2
                R_minus2=(yd+dS[j])**2+zd**2
                Rd_plus2=(y2d+dS[j])**2+z2d**2
                Rd_minus2=(y2d-dS[j])**2+z2d**2
                # Biot-Savartの法則の係数 1/2pi*h を計算 ... (8)
                Q[i,j] = -1/(2*pi)*(( (yd-dS[j])/R_plus2 - (yd+dS[j])/R_minus2 )*cos(phi[i]-phi[j])\
                    +( zd/R_plus2 - zd/R_minus2 )*sin(phi[i]-phi[j])\
                    + ( (y2d-dS[j])/Rd_minus2 - (y2d+dS[j])/Rd_plus2 )*cos(phi[i]+phi[j])\
                    +( z2d/Rd_minus2 - z2d/Rd_plus2 )*sin(phi[i]+phi[j]))
        
        # 変数を正規化
        dsigma=dS/l_e  # dS → dsigma
        nu=y/l_e  # y  → nu
        zeta=z/l_e  # z  → zeta
        q=Q*l_e  # Q  → q

        # 制約条件を定式化

        # beta : 翼根の曲げモーメントの制約条件に対するパラメータ ... (17)
        c=2*dsigma*cos(phi)  # 全揚力とスパンの荷重分布を関連付ける定数 ... (20)
        b=3*pi/2*(nu*cos(phi)+zeta*sin(phi)) * dsigma  # 翼根の曲げモーメントとスパンの荷重分布を関連付ける定数
        g=np.zeros(shape=(N))  # 揚力分布(今回最適化したい変数) ... (19)
        A=pi*q*dsigma  # ... (22)

        # これらを用いると
        # (揚力固定の条件): \bf{c}^T\cdot\bf{g}=1 ... (23)
        # (構造の制約条件): \bf{b}^T\cdot\bf{g}=beta ... (24)


        # 最適化問題 ... (26)　を解く

        # \frac{1}{e}=\bf{g}^T\cdot{A}\cdot\bf{g} ... (25)
        # が最小となる\bf{g}が求めるもの
        # Lagrangeの未定乗数法により解く

        # 係数行列 (33)式
        _A=np.c_[
            np.r_[A+A.T, np.r_[-c, -b].reshape(2, N)], np.r_[-c, 0, 0], np.r_[-b, 0, 0]]
        d=np.r_[np.zeros(shape=(N)), -1, -beta]
        # 連立N+2元一次方程式を解く
        solved_matrix=np.linalg.solve(_A, d)
        g=solved_matrix[0:N]
        # mu=solved_matrix[N+1:] # Lagrangeの未定乗数法の係数

        # 循環を求める
        # (19)式より
        # \Gamma_i=\frac{g_i}{2l_e\{rho}U}
        self.gamma=self.gamma+lr*(g*Lift/(2*l_e*rho*U)-self.gamma)

        # 吹きおろし速度を求める ... (7)
        for i in range(N):
            self.V_n[i]=np.dot(Q[i], self.gamma)

        # 揚力を計算 ... (4)
        self.local_lift=2*rho*U*self.gamma*cos(phi)*dS
        # 各点の曲げモーメントを計算
        for i in range(N):
            # i番目のパネルより右側のパネルのモーメントを考える。
            self.local_bending_moment[i]=np.dot(self.local_lift[i:],y[i:]-y[i])+np.dot(self.local_lift[i:]*tan(phi[i:]),z[i:]-z[i])

        # 誘導抗力を計算 ... (6)
        self.local_induced_drag=rho*self.gamma*self.V_n*dS
    def update(self):
        # 最適化の結果をFeed backする

        # 今処理しているパネルの上反角
        phi=0
        # 断面二次モーメント
        I=pi*((D+d)**4-D**4)/64
        y=np.linspace(self.dS[0],self.l_e-self.dS[N-1],N)
        z=np.zeros(shape=(N))
        kappa=self.local_bending_moment/(E*I)*self.dS
        for i in range(1,N):
            phi+=kappa[i-1]
            self.phi[i]=phi
            y[i]=y[i-1]+2*self.dS[i]*cos(self.phi[i])
            z[i]=z[i-1]+2*self.dS[i]*sin(self.phi[i])
        self.y=y
        self.z=z
    def plot(self):
        print("\nlift:{:2f}".format((self.local_lift).sum()*2))
        print("induced drag:{:2f}".format(self.local_induced_drag.sum()*2))
        plt.rcParams["font.size"] = 6
        _,axes=plt.subplots(3,3)
        axes[0][0].plot(self.y/cos(self.phi),self.gamma)
        axes[0][0].set_xlabel("span [m]")
        axes[0][0].set_ylabel("circulation [m^2/s]")
        axes[0][0].grid()
        axes[0][1].plot(self.y/cos(self.phi),self.local_lift/self.dS)
        axes[0][1].set_xlabel("span [m]")
        axes[0][1].set_ylabel("lift [N]")
        axes[0][1].grid()
        axes[0][2].plot(self.y/cos(self.phi),self.local_induced_drag/self.dS)
        axes[0][2].set_xlabel("span [m]")
        axes[0][2].set_ylabel("induced drag [N]")
        axes[0][2].grid()
        axes[1][0].plot(self.y/cos(self.phi),self.local_bending_moment)
        axes[1][0].set_xlabel("span [m]")
        axes[1][0].set_ylabel("bending moment [N.m]")
        axes[1][0].grid()
        axes[1][1].plot(self.y/cos(self.phi),self.V_n)
        axes[1][1].set_xlabel("span [m]")
        axes[1][1].set_ylabel("downwash velocity [m/s]")
        axes[1][1].grid()
        axes[1][2].plot(self.y/cos(self.phi),self.V_n/U)
        axes[1][2].set_xlabel("span [m]")
        axes[1][2].set_ylabel("downwash angle [rad]")
        axes[1][2].grid()
        axes[2][0].plot(self.y,self.z)
        axes[2][0].set_xlabel("y [m]")
        axes[2][0].set_ylabel("z [m]")
        axes[2][0].set_xlim(0,self.l_e)
        axes[2][0].set_ylim(0,self.l_e)
        axes[2][0].grid()
        axes[2][1].plot(self.y,self.phi)
        axes[2][1].set_xlabel("y [m]")
        axes[2][1].set_ylabel("phi [rad]")
        axes[2][1].grid()
        plt.tight_layout()
        plt.show()

wing=Wing()
for i in range(ITER):
    lr=(ITER-i)/ITER*0.1
    print("\rITER:{}/{} LOSS:{:2e}/{}".format(i+1,ITER,Lift-wing.local_lift.sum()*2,Lift),end="")
    wing.optim()
    wing.update()
wing.plot()
