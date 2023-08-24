## pythonのインストール手順

1.Windowsボタンから、スタートページを開き、右上の「すべてのアプリ」をクリックします。
2.一番下までスクロールして、「ポータルサイト」を選択します。
3.サインインしたら、「最近公開されたアプリ」の「Python」をクリックし、続けて「インストール」をクリックします。
4.python-3.11.4-amd64のインストールが完了します。


#### Pythonの動作確認
1.Windowsボタン横の検索窓から、「PowerShell」と検索し、
2.青または黒い画面が現れたら、カーソルの位置を変えずに「python」と入力して、エンターキーをおしてください。
　「Python 3.33.4 .....」という出力が現れたら、pythonが問題なくインストールされていることを意味します。
 　また、「>>>」 以降はpythonの対話型インタープリターとなっており、コードを書くことでpythonを実行できます。
  　試しに、「 print("hello") 」や「1+1」と入力して動作を確認しましょう。
3.「ctrl + z」を押した後、Enterキーを押すとPythonの入力画面から出ることが出来ます。


#### Pythonの環境構築
実際にコードを開発するためには、JupyterLabというアプリケーション(テキストエディタ)が必要です。これを環境構築と言います。 
色々な方法がありますが、下記では最も容易な環境構築の方法を紹介致します。

1.Windowsボタン横の検索窓から、「PowerShell」と検索し、エンターキーを押します。
2.青または黒い画面が現れたら、カーソルの位置を変えずに下記コードを実行しましょう。

~~~
python -m venv env
~~~

~~~
cd env\Scripts
~~~

~~~
pip install jupyterlab --trusted-host pypi.python.org --trusted-host files.pythonhosted.org --trusted-host pypi.org
~~~


実行すると、JupyterLabのインストールが開始します。完了するまで、暫く待機します。


4.インストール完了したら、
~~~
WARNING: The script ***.exe  is installed in 'C:\Users\H~~~~\AppData\Roaming\Python\Python311\Scripts
~~~
という出力の　''　で囲まれた部分をコピーし、検索窓にペーストして、そのパスを開きます。  
そこに、jupyter-lab.exeというファイルがあれば、そのファイル上で右クリックし、「ショートカットの作成」を行います。
デスクトップ上にショートカットを作成し、そこからjupyter-labを実行できるようにします。


#### JupyterLabの起動

jupyter-labを実行し、起動します。