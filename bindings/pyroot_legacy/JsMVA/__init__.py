# -*- coding:utf-8 -*-
## @mainpage
# @package JsMVA
# @author  Attila Bagoly <battila93@gmail.com>

from IPython.core.extensions import ExtensionManager

## This function will register JsMVAMagic class to ipython
def loadExtensions():
    ip     = get_ipython()
    extMgr = ExtensionManager(ip)
    extMgr.load_extension("JsMVA.JsMVAMagic")

loadExtensions()