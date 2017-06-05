import socket
import unittest


fqdn = socket.getfqdn()
if not fqdn[-6:] == 'uib.no':
	raise unittest.SkipTest('Tests rely on UiB infrastructure')

# C'est le fin
