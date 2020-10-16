#include "UpdateBroadcaster.hpp"

void UpdateBroadcaster::subscribe(IPredictor *subscriber) {
  subscribers[n] = subscriber;
  n++;
  assert(n < 1024);
}

void UpdateBroadcaster::broadcastUpdate() { //广播更新
  for( int i = 0; i < n; i++ ) {
    IPredictor *subscriber = subscribers[i];
    subscriber->update();
  }
  n = 0;
}
