D=load('tabledata05.txt');
loglog(D([3:9],1),D([3:9],2));hold on
(log(D(9,2))-log(D(3,2)))/(log(D(9,1))-log(D(3,1)))

loglog(D([3:9],1),D([3:9],3));
(log(D(9,3))-log(D(3,3)))/(log(D(9,1))-log(D(3,1)))

loglog(D([3:9],1),D([3:9],4));
(log(D(9,4))-log(D(3,4)))/(log(D(9,1))-log(D(3,1)))