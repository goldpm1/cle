import scipy
import numpy as np

def main(mixture):
    NUM_BLOCK = mixture.shape[0]
    zeropoint = np.zeros(NUM_BLOCK)

    new_mixture  =mixture 

    while True:
        new_mixture =  (new_mixture.transpose()[new_mixture.transpose()[:, 0].argsort()].transpose())       # block 0의 mixutre 크기를 기준으로 오름차순 정렬

        min_distance = float("inf")
        for i in new_mixture.transpose()[1:,]:           # 정렬해서 맨 앞 clone을 빼고 나머지를 돌고, min_distance를 구함
            min_distance = np.min([min_distance, scipy.spatial.distance.euclidean(zeropoint, i)])
        
        d0 = scipy.spatial.distance.euclidean(zeropoint, new_mixture[:, 0])

        if d0 < min_distance:      # 맨 앞 clone이 거리적으로 나머지 모든 clone보다 안쪽이라고 생각될 때
            new_mixture = new_mixture[:, 1:]      # 맨 앞 clone을 빼준다
        else:
            break                  # 더 이상 그런 관계가 없으면 종료해준다

    return new_mixture