# kaggle
My Kaggle projects.

## Titanic

This is the [Titanic - Machine Learning from Disaster](https://www.kaggle.com/competitions/titanic) competition.

Afted data cleaning ([titanic_data_cleaning.ipynb](https://github.com/tmaciazek/kaggle/blob/main/titanic/titanic_data_cleaning.ipynb)), I approach this using three techniques:

-  logistic regression using polynomial features (of order $2$) combined with PCA ([logistic.ipynb](https://github.com/tmaciazek/kaggle/blob/main/titanic/logistic.ipynb)),
-  a single decision tree ([tree.ipynb](https://github.com/tmaciazek/kaggle/blob/main/titanic/tree.ipynb)),
-  a random forest ([random_forest.ipynb](https://github.com/tmaciazek/kaggle/blob/main/titanic/random_forest.ipynb)).

## Natural Language Processing with Disaster Tweets

This is the [Natural Language Processing with Disaster Tweets](https://www.kaggle.com/competitions/nlp-getting-started) competition.

The file [tweets_preprocessing.ipynb](https://github.com/tmaciazek/kaggle/blob/main/NLP_disaster_tweets/tweets_preprocessing.ipynb) contains data cleaning and processing, in particular:
-  inferring location country from location data and from the text of the tweet,
-  tweet text cleaning and tokenization,
-  applying word embeddings using the GloVe dataset [glove.twitter.27B](https://nlp.stanford.edu/projects/glove/).

The cleaned data is subsequently fed into a recurrent neural net with two bidirectional LSTM layers shown in the schamtic picture below.

<img src="NLP_disaster_tweets/net.png" width="620" height="500"> 

This can be found in the notebook [lstm.ipynb](NLP_disaster_tweets/lstm.ipynb).
