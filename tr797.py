import numpy as np
from numpy import cos,sin,pi,sqrt
from matplotlib import pyplot as plt

# 定数（既知とするパラメータ）
N = 100 # 分割数
l_e = 15.0 # 右翼の長さ
dS = np.ones(shape=(N))*l_e/(N*2) # パネルの半幅
rho = 1.2  # 空気の密度 [kg/m^3]
U = 7.2  # 一様流の流速 [m/s]
Lift=1000.0 # 揚力 [N]
beta=0.85


# Trefftz面での翼の形
# TODO 曲げモーメントによって平面形を変化させる
y = np.linspace(dS[0], l_e-dS[N-1], N)
z = np.zeros(shape=(N))
# 上反角の分布
phi = np.zeros(shape=(N))


# Biot-Savartの法則の係数h/2pi
Q = np.zeros(shape=(N, N))
q = np.zeros(shape=(N, N))

# 翼の形状の定式化
for i in range(N):
    for j in range(N):
        # (9)式 ... 回転した座標系における各パネルの座標
        # 右側の翼について
        ydash = (y[i]-y[j])*cos(phi[j])+(z[i]-z[j])*sin(phi[j])
        zdash = -(y[i]-y[j])*sin(phi[j])+(z[i]-z[j])*cos(phi[j])
        # 左側の翼について( y → -y , phi[j] → -phi[j])
        y2dash = (y[i]+y[j])*cos(phi[j])-(z[i]-z[j])*sin(phi[j])
        z2dash = (y[i]+y[j])*sin(phi[j])+(z[i]-z[j])*cos(phi[j])

        # (10)式 ... 注目する２つのパネルの距離（plus→左側、minus→右側）
        # 右側の翼について
        R_plus2 = (ydash-dS[j])**2+zdash**2
        R_minus2 = (ydash+dS[j])**2+zdash**2
        # 左側の翼について
        Rdash_plus2 = (y2dash+dS[j])**2+z2dash**2
        Rdash_minus2 = (y2dash-dS[j])**2+z2dash**2

        # (8)式
        # 循環により誘導される速度の係数 h/2piを計算
        Q[i,j] = -1/(2*pi)*(( (ydash-dS[j])/R_plus2 - (ydash+dS[j])/R_minus2 )*cos(phi[i]-phi[j])\
                +( zdash/R_plus2 - zdash/R_minus2 )*sin(phi[i]-phi[j])\
                + ( (y2dash-dS[j])/Rdash_minus2 - (y2dash+dS[j])/Rdash_plus2 )*cos(phi[i]+phi[j])\
                +( z2dash/Rdash_minus2 - z2dash/Rdash_plus2 )*sin(phi[i]+phi[j]))

# 変数を正規化

dsigma=dS/l_e  # dS → dsigma
nu=y/l_e  # y  → nu
zeta=z/l_e  # z  → zeta
q=Q*l_e  # Q  → q

# 制約条件を定式化

# beta : 翼根の曲げモーメントの制約条件に対するパラメータ ... (17)
c=2*dsigma*cos(phi)  # 全揚力とスパンの荷重分布を関連付ける定数 ... (20)
b=3*pi/2*(nu*cos(phi)+zeta*sin(phi)) * \
             dsigma  # 翼根の曲げモーメントとスパンの荷重分布を関連付ける定数
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
mu=solved_matrix[N+1:]

# 循環を求める
# (19)式より
# \Gamma_i=\frac{g_i}{2l_e\{rho}U}
Gamma=g*Lift/(2*l_e*rho*U)

# 吹きおろし速度を求める ... (7)
V_n=np.zeros(shape=(N))
for i in range(N):
    V_n[i]=np.dot(Q[i], Gamma)

# 単位長さ当たりの揚力を計算 ... (4)
local_lift=2*rho*U*Gamma*cos(phi)*dS
# 各点の曲げモーメントを計算
local_bending_moment=np.zeros(shape=(N))
for i in range(N):
    # i番目のパネルより右側のパネルのモーメントを考える。
    tmp_b_mom=0
    for j in range(i,N):
        tmp_b_mom+=local_lift[j]*(y[j]-y[i])
    local_bending_moment[i]=tmp_b_mom
# 単位長さ当たりの誘導抗力を計算 ... (6)
local_induced_drag=2*rho*Gamma*V_n*dS


# グラフを表示
plt.rcParams["font.size"] = 6
fig,axes=plt.subplots(2,3)
axes[0][0].plot(y,Gamma)
axes[0][0].set_xlabel("span[m]")
axes[0][0].set_ylabel("circulation")
axes[0][0].grid()
axes[0][1].plot(y,local_lift)
axes[0][1].set_xlabel("span[m]")
axes[0][1].set_ylabel("lift[N]")
axes[0][1].grid()
axes[0][2].plot(y,local_induced_drag)
axes[0][2].set_xlabel("span[m]")
axes[0][2].set_ylabel("induced drag[N]")
axes[0][2].grid()
axes[1][0].plot(y,local_bending_moment)
axes[1][0].set_xlabel("span[m]")
axes[1][0].set_ylabel("bending moment[Nm]")
axes[1][0].grid()
axes[1][1].plot(y,V_n)
axes[1][1].set_xlabel("span[m]")
axes[1][1].set_ylabel("induced vertial velocity[m/s]")
axes[1][1].grid()
# axes[2][0].set_title("Lift")
# axes[2][0].plot(y,L)
# axes[2][1].set_title("Lift")
# axes[2][1].plot(y,L)
# axes[2][2].set_title("Lift")
# axes[2][2].plot(y,L)
plt.show()