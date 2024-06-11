from sklearn.datasets import make_circles
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
import numpy as np

# 生成一个带有噪声的二维圆形数据集
X, y = make_circles(n_samples=200, noise=0.1, factor=0.5, random_state=42)

# 可视化数据集
plt.scatter(X[:,0], X[:,1], c=y, cmap='bwr', alpha=0.5)
plt.title("Generated Circular Dataset")
plt.show()

# 定义三种核函数
kernels = ['linear', 'poly', 'rbf']

# 分别用三种核函数训练SVM分类器
for kernel in kernels:
    clf = SVC(kernel=kernel)
    clf.fit(X, y)
    y_pred = clf.predict(X)
    accuracy = accuracy_score(y, y_pred)
    print(f"Accuracy with {kernel} kernel: {accuracy:.2f}")
    
    # 可视化决策边界
    h = 0.02
    x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
    y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    plt.contourf(xx, yy, Z, cmap=plt.cm.coolwarm, alpha=0.8)
    plt.scatter(X[:,0], X[:,1], c=y, cmap='bwr', edgecolors='k')
    plt.title(f"Decision Boundary with {kernel} Kernel")
    plt.show()