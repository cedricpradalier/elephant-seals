
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.neural_network import MLPRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score

SexMap={'M':0, 'F':1}
PoseMap={'pose':0, 'F':1}

try:
    plt.close()
except:
    pass

def scaling_da(X,y):
    scale=1.+np.random.randn()*0.05
    return np.hstack([X[:-2]*scale,X[-2:]]),y*(scale**3)

def generate_augmented_ds(X,y,n_new):
    Ox=[X]
    if len(y.shape)==1:
        y=np.array([y]).T
    Oy=[y]
    for i in range(n_new):
        irand=np.random.randint(X.shape[0])
        Xnew,ynew=scaling_da(X[irand,:],y[irand,0])
        Ox.append(np.array([Xnew]))
        Oy.append(np.array([ynew]))
    return np.vstack(Ox),np.vstack(Oy)

        

df = pd.read_csv("mensurations_elephs.csv")
df['sexe_num'] = df['sexe'].apply(lambda s: SexMap[s])
df['pose_num'] = df['pose'].apply(lambda s: 0 if s=='pose' else 1)

input_col=['longueur', 'circonfa', 'circonfb', 'circonfc', 'circonfd', 'distance_a_museau',
        'distance_a_b', 'distance_b_c', 'distance_c_d', 'distance_d_queue', 'sexe_num', 'pose_num']

df_learn=df.loc[:,['masse']+input_col]
df_learn.dropna(axis=0, how='any', inplace=True)
y=df_learn['masse'].values
ym=y.mean()
ystd=y.std()
y=(y-ym)/ystd

X=df_learn[input_col].values
Xm=X[:,:10].mean(axis=0)
Xstd=X[:,:10].std(axis=0)
X[:,:10]=(X[:,:10]-Xm)/Xstd
X=np.hstack([X,np.ones((X.shape[0],1))])

df_all=df.dropna(axis=0, how='any', inplace=False, subset=input_col).copy(deep=True)

X_all=df_all.loc[:,input_col].values
X_all[:,:10]=(X_all[:,:10]-Xm)/Xstd
X_all=np.hstack([X_all,np.ones((X_all.shape[0],1))])

# df_all.loc[:,'masse_predite']=y_all
# df_all.to_csv("mensurations_elephs_pred.csv",na_rep='NULL',index=False,float_format="%02g")

bpdata=[]
bplabel=[]
bpmodels=[]

def pscore(y_ref, y_pred, label):
    global ym,ystd
    y_ref = y_ref*ystd+ym
    y_pred = y_pred*ystd+ym
    score=np.sqrt(np.sum((y_ref-y_pred)**2)/y_ref.shape[0])
    print("RMSE %f" % score)
    score=(np.sum(np.abs(y_ref-y_pred))/y_ref.shape[0])
    print("MAE %f" % score)
    percentage_error=100*(y_ref-y_pred)/y_ref
    score=(np.sum(np.abs(percentage_error))/y_ref.shape[0])
    print("MAPE %f" % score)
    bpdata.append(percentage_error)
    bplabel.append(label)


def fit():
    global bpmodels
    bpmodels=[]
    if False:
        Xaug,yaug=generate_augmented_ds(X,y,10000)
        X_train, X_test, y_train, y_test = train_test_split(Xaug, yaug.flatten(), random_state=1, test_size=0.25)
    else:
        X_train_org, X_test, y_train_org, y_test = train_test_split(X, y, random_state=1, test_size=0.1)
        X_train,y_train=generate_augmented_ds(X_train_org,y_train_org,10000)
        y_train=y_train.flatten()

    # Linear regression
    print("Linear model All")
    A=np.hstack([X_train,np.ones((X_train.shape[0],1))])
    B=np.array([y_train]).T
    model=np.linalg.pinv(A) @ B
    bpmodels.append(model)
    y_pred0=np.hstack([X_test,np.ones((X_test.shape[0],1))]) @ model
    y_pred0=y_pred0.flatten()
    pscore(y_test,y_pred0,"Linear all")

    print("Linear model Female")
    idxf=np.where(X_train[:,10]==1)[0]
    A=np.hstack([X_train[idxf,:],np.ones((idxf.shape[0],1))])
    B=np.array([y_train[idxf]]).T
    model_female=np.linalg.pinv(A) @ B
    bpmodels.append(model_female)
    idxf=np.where(X_test[:,10]==1)[0]
    y_pred0f=np.hstack([X_test[idxf,:],np.ones((idxf.shape[0],1))]) @ model_female
    y_pred0f=y_pred0f.flatten()
    pscore(y_test[idxf],y_pred0f,"Linear female")

    print("Linear model Male")
    idxm=np.where(X_train[:,10]==0)[0]
    A=np.hstack([X_train[idxm,:],np.ones((idxm.shape[0],1))])
    B=np.array([y_train[idxm]]).T
    model_male=np.linalg.pinv(A) @ B
    bpmodels.append(model_male)
    idxm=np.where(X_test[:,10]==1)[0]
    y_pred0m=np.hstack([X_test[idxm,:],np.ones((idxm.shape[0],1))]) @ model_male
    y_pred0m=y_pred0m.flatten()
    pscore(y_test[idxm],y_pred0m,"Linear male")

    print("Linear model Combined")
    pscore(np.array(list(y_test[idxm])+list(y_test[idxf])),
            np.array(list(y_pred0m)+list(y_pred0f)),"Linear combined")
    bpmodels.append((model_female,model_male))

    print("MLP model")
    mlp_regr = MLPRegressor(random_state=1, max_iter=2000, tol=0.1, hidden_layer_sizes=(64,64,),learning_rate='adaptive',activation='relu')
    mlp_regr.fit(X_train, y_train)
    bpmodels.append(mlp_regr)
    y_pred1=mlp_regr.predict(X_test)
    print(mlp_regr.predict(X_test[:2]))
    print(mlp_regr.score(X_test, y_test))
    pscore(y_test,y_pred1,"MLP")
    # score=cross_val_score(regr, X, y, cv=10, n_jobs=10)
    # print(score)

    print("DecisionTree model")
    dt_regr = DecisionTreeRegressor(random_state=1,max_depth=32,criterion='mse')
    dt_regr.fit(X_train, y_train)
    bpmodels.append(dt_regr)
    y_pred2=dt_regr.predict(X_test)
    print(dt_regr.predict(X_test[:2]))
    print(dt_regr.score(X_test, y_test))
    pscore(y_test,y_pred2,"Decision Tree")
    # score=cross_val_score(regressor, X, y, cv=10, n_jobs=10)
    # print(score)

    print("RandomForest model")
    rf_regr = RandomForestRegressor(max_depth=32, random_state=1, criterion='mse', n_jobs=16)
    rf_regr.fit(X_train, y_train)
    bpmodels.append(rf_regr)
    y_pred3=rf_regr.predict(X_test)
    print(rf_regr.predict(X_test[:2]))
    print(rf_regr.score(X_test, y_test))
    pscore(y_test,y_pred3,"Random Forrest")

    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.set(axisbelow=True,  # Hide the grid behind plot objects
            title='Comparison of mass prediction models',
            xlabel='Models',
            ylabel='Prediction error in percent',
    )
    ax1.boxplot(bpdata)
    ax1.set_xticklabels(bplabel, rotation=45, fontsize=8, verticalalignment='center')
    ax1.grid()
    plt.savefig("models.png")

    y_all = mlp_regr.predict(X_all)*ystd+ym

    pred=np.array([y_test,y_pred0,y_pred1,y_pred2,y_pred3]).T
    pred=pred * ystd + ym
    err=100* (pred - np.array([pred[:,0]]).T @ np.array([[1]*pred.shape[1]]))/pred
    print(pred)
    print("RMSE: ")
    print(np.sqrt(np.power(err,2).mean(axis=0)))
    print("MAE: ")
    print(np.abs(err).mean(axis=0))
    return pred,err,y_all

