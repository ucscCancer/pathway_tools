#!/usr/bin/env python

import tornado.ioloop
import tornado.web
import tornado.httpserver
import urllib
import datetime
import base64
import time
import json
import re
import time
import sys
from copy import deepcopy

if __name__ == "__main__":
    application = tornado.web.Application([
        (r"/(.*)", tornado.web.StaticFileHandler, {"path": sys.argv[1]}),
    ])
    http_server = tornado.httpserver.HTTPServer(application)
    http_server.listen(sys.argv[2])
    tornado.ioloop.IOLoop.instance().start()
