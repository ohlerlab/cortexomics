
if __name__ = '__main__':
	def add1(x):
		return(x+1)

	pool = multiprocessing.Pool(10)


	foo = 1
	fooa = pool.apply_async(add1,args=(foo,))
	fooa.get()

	pool.map(add1, range(10))

	#create a list of threads
	urls = ['a','b','c']
	def crawl(s,result,index):
		
	threads = []
	results = [{} for x in urls]
	# In this case 'urls' is a list of urls to be crawled.
	for ii in range(len(urls)):
	    # We start one thread per url present.
	    process = Thread(target=crawl, args=[urls[ii], result, ii])
	    process.start()
	    threads.append(process)
	# We now pause execution on the main thread by 'joining' all of our started threads.
	# This ensures that each has finished processing the urls.
	for process in threads:
	    process.join()




	from concurrent.futures import ProcessPoolExecutor

	num_processes = 3

	def add1(x):
		return(x+1)


	list(results)
	import multiprocessing

	multiprocessing.cpu_count()

############