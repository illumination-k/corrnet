# TODO List

- [x] [conetk/models/expression/network.py](https://github.com/sepro/CoNekT/blob/master/conekt/models/expression/networks.py)の中のlrstripを読む関数の中身を使ってhighest reciprocal rankのテストをしておく。~~なんか逆側を見ていて、逆側が閾値以下ならNoneにして、Noneとrankがあって、rankが閾値以下ならrankをとってきているっぽいな。~~多分嘘。少なくとも双方向の情報を両方保持しているので、数は二倍程度になるっぽい。あとはTranscriptsの差があるので、たぶんOK!性染色体のデータどうするかかな。
- [] rankは全部f64でいい
- [] serializeを`&[u8]`でとったほうがよいのでそのあたりを全部修正する
- [] hcca.pyのrust実装
- [] cluster間のjaccad index
- [x] codon usageの計算
- [] logit scoreの実装
    - [] 論文読む
    - [] 実装する
