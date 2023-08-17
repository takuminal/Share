"""Erastic Net による機械学習と、
ベイズ最適化によるパラメータチューニング、
交差検証による汎化性能確認をまとめて関数にした。
"""
def Elastic_Net(properties,X,y,test_train_ratio,N_trial):

    from sklearn.linear_model import ElasticNet
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    import optuna
    import logging
    #出力のログを制限
    optuna.logging.set_verbosity(optuna.logging.WARNING)


    from sklearn.model_selection import train_test_split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_train_ratio,random_state=42)

    
    from sklearn.model_selection import KFold
    from sklearn.model_selection import cross_val_score

    alpha =0.1
    rho =0.5
    
    
    kfold = KFold(n_splits=5, shuffle=True, random_state=42) # 交差検証のためのk-foldオブジェクトを作成
    model = ElasticNet(alpha=alpha, l1_ratio=rho, max_iter=10000)
    scores = cross_val_score(model, X, y, cv=kfold) # k-分割交差検証を実行し、モデルの性能を評価


    #エラスティックネットのハイパラ->スコアの関数を作成。
    def objective(trial):
        #ハイパーパラメータの設定
        alpha = trial.suggest_float("alpha",0,1)
        rho = trial.suggest_float("rho",0,1)
        kfold = KFold(n_splits=3, shuffle=True, random_state=42) # 交差検証のためのk-foldオブジェクトを
        model = ElasticNet(alpha=alpha, l1_ratio=rho, max_iter=10000)
        scores = cross_val_score(model, X, y, cv=kfold) # k-分割交差検証を実行し、モデルの性能を評価

         # 平均評価スコアを計算して表示
        return scores.mean()

    #optunaでの最適化
    study = optuna.create_study(direction='maximize')
    study.optimize(objective, n_trials = N_trial)
    print(f"best value of parameters : {study.best_params}")


    #最適化パラメータでの再推論
    alpha = study.best_params["alpha"]  # 正則化の強度 (alpha > 0)
    rho = study.best_params["rho"]  # Elastic Netのハイパーパラメータ (0 <= rho <= 1)
    elastic_net = ElasticNet(alpha=alpha, l1_ratio=rho, max_iter=10000)
    elastic_net.fit(X_train, y_train)
    y_pred = elastic_net.predict(X_test)
    en_score = elastic_net.score(X_test,y_test)
    print(f"Test Score via Optimized Parameters : {en_score}")

    #EN回帰の係数
    EN_coef = elastic_net.coef_

    # 試行の履歴がdictのリストで得られる。
    n = [i.number for i in study.trials]
    score_ls = [i.values for i in study.trials]

 
    #結果の可視化
    plt.figure(figsize=(12,8))
    plt.style.use("_mpl-gallery-nogrid")


    plt.subplot(2, 1, 1)  # (行数, 列数, サブプロット番号)
    #軸のフォントの調整
    plt.xticks([])
    plt.yticks(fontsize=10)

    #絶対値が最大の要素の1.5倍の正負がyの表示領域
    abs_max=np.max(np.abs(EN_coef))
    plt.ylim(-abs_max*1.5,abs_max*1.5)

    plt.title("Model Coefficient")
    plt.xlabel("wavenumber(cm-1)")
    plt.ylabel("coefficient")
    plt.grid(True)

    for i, wn in enumerate(properties):
        if EN_coef[i] > 0:
            plt.text(wn,EN_coef[i],str(wn),
                     rotation=90,
                     ha='center',
                     va='bottom')
        elif EN_coef[i] < 0:
            plt.text(wn,EN_coef[i],str(wn),
                     rotation=90,
                     ha='center',
                     va='top')

    plt.bar(properties,EN_coef,color ="lightblue")


    plt.subplot(2, 2, 3)
    # サブプロット1のプロット
    plt.scatter(y_test,y_pred,color="gray")
    # 直線 y = x を点線でプロット
    plt.plot(y_test, y_test, linestyle='dashed', color='aquamarine')
    plt.xlabel("y_test")
    plt.ylabel("y_predict")
    plt.title("ElasticNet CrossVaridation")
    plt.grid(True)
    # テキストを表示
    text_x = np.min(y_test)  # テキストのx座標
    text_y = np.max(y_test)  # テキストのy座標


    text = f"Best Score: {elastic_net.score(X_test,y_test).round(2)}" # 表示するテキスト
    plt.text(text_x, text_y, text, fontsize=12, ha='left', va='top', color='k')

    # 2行目に2列のサブプロットを作成（右側）
    plt.subplot(2, 2, 4)  # (行数, 列数, サブプロット番号)
    # サブプロット3のプロット
    plt.plot(n,score_ls,color ="gray")
    # 直線 y = 1 を点線でプロット
    plt.fill_between(n,np.ones(len(n)),np.ones(len(n))*0.8,color="aquamarine",alpha=0.7)
    plt.xlabel("bayes optimization n")
    plt.ylabel("test score")
    plt.title("Hyper parameter Optimize")
    plt.grid(True)

    plt.show()

