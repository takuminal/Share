import subprocess

# パッケージのリストを指定してpip installを実行する
package_list = ['pubchempy',
                'PySimpleGUI',
                'rdkit',
                "sklearn",
                "matplotlib",
                "numpy",
                "optuna",
                "logging"]
subprocess.run(['pip', 'install'] + package_list)

import pkg_resources

def is_package_installed(package_name):
    try:
        pkg_resources.get_distribution(package_name)
        return True
    except pkg_resources.DistributionNotFound:
        return False

# パッケージ名を指定して確認する
checker = True

for package_name in  package_list:

    if is_package_installed(package_name):
        print(f"{package_name} finished")
        checker = checker*True
    else:
        print(f"{package_name} error")
        checker = checker*False

if checker ==True:
    print("■■        ■■        ■")
    print(" ■   ■■   ■         ■")
    print(" ■   ■■   ■   ■■■   ■    ■■■■  ■■■■   ■■■■■■■■    ■■■")
    print(" ■   ■■   ■  ■■  ■  ■   ■■    ■■  ■   ■   ■   ■  ■■  ■")
    print(" ■■ ■■ ■  ■  ■   ■■ ■   ■     ■    ■  ■   ■   ■  ■   ■■")
    print(" ■■ ■  ■ ■■  ■■■■■■ ■   ■     ■    ■  ■   ■   ■  ■■■■■■")
    print("  ■ ■  ■ ■   ■      ■   ■     ■    ■  ■   ■   ■  ■     ")
    print("  ■■    ■■   ■■     ■   ■■    ■■  ■   ■   ■   ■  ■■    ")
    print("  ■■    ■■    ■■■■   ■   ■■■■  ■■■■   ■   ■   ■   ■■■■ ")
import time
# 5秒間待機
time.sleep(10)