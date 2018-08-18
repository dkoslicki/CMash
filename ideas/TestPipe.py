import multiprocessing

def worker(procnum, send_end):
	'''worker function'''
	for i in range(10):
		result = str(i) + ' represent!'
		#print result
		send_end.send(result)

def main():
	jobs = []
	pipe_list = []
	for i in range(5):
		recv_end, send_end = multiprocessing.Pipe(False)
		p = multiprocessing.Process(target=worker, args=(i, send_end))
		jobs.append(p)
		pipe_list.append(recv_end)
		p.start()

	for proc in jobs:
		proc.join()
	result_list = [x.recv() for x in pipe_list]
	print result_list

if __name__ == '__main__':
	main()