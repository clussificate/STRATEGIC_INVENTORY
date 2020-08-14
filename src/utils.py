# -*- coding: utf-8 -*-
"""
@Created at 2020/8/14 16:56
@Author: Kurt
@file:utils.py
@Desc:
"""


class DefinedException(Exception):
    def __init__(self, err=None):
        self.err = err
        Exception.__init__(self, err)