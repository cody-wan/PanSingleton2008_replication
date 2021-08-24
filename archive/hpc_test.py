import time
from datetime import datetime
import multiprocessing

def func(i):
    with open("test.txt", "a") as f:
        f.write(f"{i} - {datetime.now()}\n")
    time.sleep(5)
    with open("test.txt", "a") as f:
        f.write(f"{i} - {datetime.now()}\n")


if __name__ == "__main__":
    ps = [multiprocessing.Process(target=func, args=(i,)) for i in range(60)]
    for p in ps:
        p.start()
    for p in ps:
        p.join()
