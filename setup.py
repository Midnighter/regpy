#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
===================================
RegPy - Holistic Genetic Regulation
===================================

:Authors:
    Moritz Emanuel Beber
:Date:
    2011-11-23
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    setup.py
"""


from distutils.core import setup


setup(
    name = "regpy",
    version = "0.1",
    description = "classes and functions that model genetic regulation",
    author = "Moritz Emanuel Beber",
    author_email = "moritz (dot) beber (at) googlemail (dot) com",
    url = "http://github.com/Midnighter/",
    packages = [
            "regpy",
            "regpy.model",
            ],
    )

