import time


start  = time.time()

# for i in range(10000000):
#     a = False
#     if 1 > 2:
#         a = True
#     elif 2 < 3:
#         a = True
# print(time.time() - start)
#
# start  = time.time()
# for i in range(10000000):
#     if (1 > 2 or 2 < 3):
#         pass
# print(time.time() - start)


b = 3
a = 2

start  = time.time()
for i in range(10000000):
    if (b >= a):
        a = b
print(time.time() - start)


start  = time.time()
for i in range(10000000):
    a = min(b,a)
print(time.time() - start)




