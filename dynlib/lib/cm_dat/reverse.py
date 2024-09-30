#!/usr/bin/env python
# -*- encoding: utf-8

def reverse_cm(cm_data):
	cm_r = {}
	for channel, points in cm_data.items():
		cm_r[channel] = []
		for point in points[::-1]:
			cm_r[channel].append((1.0-point[0], )+point[1:])

	return cm_r

#
