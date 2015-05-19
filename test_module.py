# -*- coding: utf-8 -*-
"""
Created on Wed May 13 16:28:55 2015

@author: Hassan
"""

if __name__ == '__main__':
    print 'main'
    
A = 3

def greet():
    print "hello"
    
greet()

def f1(string):
    print(string.upper())
    
def superf():    
    def f1(string):
        print(string.lower())
        
    f1('Bazinga!')
    
superf()