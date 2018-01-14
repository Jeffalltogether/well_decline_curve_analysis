# grad descent test
import numpy as np

def decline_curve(t, qi, b, di):
	# unit of t is in days
	print qi, b, di
	return qi*(1.0 - b*(di)*t)**(-1.0 / b)

def decline_curve_residuals(startParams, t, y):
	# unit of t is in days
	q = startParams
	return q[0]*(1.0 - q[1]*(q[2])*t)**(-1.0 / q[1])

def gradient_descent(t, y, startParams, epochs = 100, learningRate = 0.01):
	import numpy as np
	mean_t = np.mean(t)
	mean_y = np.mean(y)
	t = t - mean_t
	y = y - mean_y

	qi, b, di = startParams[0], startParams[1], startParams[2]
	for epoch in range(epochs):
		# Compute predictions with parameters
		print qi
		print b
		print di

		y_pred = decline_curve(t, qi, b, di)

		# compute error
		SSE = 0.5* np.sum(y - y_pred)**2
		print 'error based on parameters %.0f, %.4f, %.4f = %.4f' %(qi, b, di, SSE)
		# Calculate gradients
		dSSE_dqi = np.sum(-1.0*(-qi*(-b*di*y + 1.0)**(-1.0/b) + t)*(-b*di*y + 1.0)**(-1.0/b))

		dSSE_db = np.sum(-1.0*qi*(-qi*(-b*di*y + 1.0)**(-1.0/b) + t)*(-b*di*y + 1.0)**(-1.0/b)*(1.0*di*y/(b*(-b*di*y + 1.0)) + 1.0*np.log10(-b*di*y + 1.0)/b**2))

		dSSE_ddi = np.sum(-1.0*qi*y*(-qi*(-b*di*y + 1.0)**(-1.0/b) + t)*(-b*di*y + 1.0)**(-1.0/b)/(-b*di*y + 1.0))
		
		# update parameters
		print dSSE_dqi
		print dSSE_db
		print dSSE_ddi

		qi = qi - learningRate * np.real(dSSE_dqi)
		b = b - learningRate * np.real(dSSE_db)
		di = di - learningRate * np.real(dSSE_ddi)


	return qi, b, di

if __name__ == '__main__':

	startParams = np.array([5000.0, 0.04, 0.])
	qi, b, di = startParams[0], startParams[1], startParams[2]


	a = np.arange(0, 510, 10)
	b = decline_curve(a, qi, b, di)
	mean_b = np.mean(b)
	st_b = np.std(b)
	# c = b + np.random.normal(mean_b, st_b)
	
	print a
	print b
