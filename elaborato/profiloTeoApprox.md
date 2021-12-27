# Esercizio teorie approssimate

## 1. Profilo F-104

Il profilo normalizzato dello *starfighter* è descritto dalla seguente curva cartesiana:  
$$
\begin{align}
y&=\pm \big[- \big(R - \frac{\tau}{2}\big) + \sqrt{R^2 - (x-0.5)^2}\big]\\
R&=\big(\frac{1}{\tau} + \tau \big)\frac c 4 \simeq 1.7 \\
\tau&=0.15\\
c&=1
\end{align}
$$  

Data la geometria del profilo è possibile applicare le seguenti teorie approssimate, dato l' elevato numero di Mach e la piccola curvatura del profilo:

* Newton
* Newton-Bousemman
* cono-tangente 
* urto-espansione

Nel codice la geometria del profilo è discretizzata in un numero finito di punti $(x,y)$. In tal modo é possibile approssimare il profilo con una linea spezzata chiusa ed ogni tratto $i-esimo$  é caratterizzato dalle seguenti caratteristiche geometriche:

- $P_i=(x_i,y_i)$ 
- $l_i= \sqrt{(x_{i+1} - x_i)^2 + (y_{i+1} - y_i)^2}$
- $\theta_{i_{geom}}= \arctan(\frac{y_{i+1} - y_i}{x_{i+1} - x_i})$

Conviene infittire il numero di punti in prossimità dei bordo di attacco e di uscita dove è maggiore la $\frac{dy} {dx}$ . Nel codice matlab è stata usata una distribuzione di punti alla Cebicev implementata nel seguente modo

```matlab
angle = linspace(pi,0,100);
x=(1+cos(angle))/2
```  

Sul bordo d'attacco l'angolo $\theta_{geom}$  sul dorso calcolato analiticamente vale:
$$
\theta_{geom}= \arctan( \frac{0.5}{\sqrt{R^2 - 0.25}}) \simeq 17°
$$
usando uno stencil di punti come prima citato il primo tratto della spezzata sul bordo d'attacco ha un angolo $\theta_{1_{geom}} \simeq 17°$

### Metodo di Newton

Secondo tale metodo il coefficiente di pressione del profilo immerso nel "fluido" é dato dalla relazione:
$$
c_p = 2\sin^2(\theta)
$$
dove $\theta$  é l'angolo che forma il punto del profilo con la corrente asintotica. quindi per un profilo posto ad un certo angolo d'attacco $\alpha$  si ha $c_p=c_p = 2\sin(\theta_{geom} \pm \alpha)^2$ . 
Dato che in questo modello la pressione é data dallo scambio di quantità di moto delle molecole con il corpo,  per le zone in cui si ha $\theta_{geom} \pm \alpha < 0$ le molecole non impattano sul corpo e quindi si creano delle zone d'ombra aerodinamica.

L' implementazione per il calcolo del cp con questo metodo é piuttosto semplice.  infatti basta semplicemente applicare la formula del $c_p$ al vettore dei $\theta_i$ prima del dorso e poi del ventre stando attenti a considerare la zona d'ombra. Per esempio per il dorso:

```matlab
cpDorsoN=2*sin( thetaDorso(thetaDorso>0)).^2;
```

### Metodo di Newton-Buseman

per tener conto della curvatura del profilo nel calcolo del coefficiente di pressione con il metodo di Newton, Buseman ha proposto la correzione di eq. (1) , dove vale il segno $+$ per corpi concavi e $-$ per corpi convessi.
$$
c_{p_B}=c_{p_N} \pm \frac 2 R \int_0^y \cos(\theta) dy \tag{1}
$$
per $\theta << 1$ ,  cioè per corpi sottili con curvatura non grande, $\cos(\theta)\simeq 1$ e la formula di prima può essere approssimata come:
$$
c_{p_B}=c_{p_N} \pm \frac 2 Ry
$$
Inoltre dato che le ipotesi dell' interazione fluido struttura sono le stesse di quelle della teoria di Newton, non hanno senso valori di $c_p < 0$. 

### Metodo del cuneo tangente

Secondo tale teoria il $c_p$ nel generico punto  $(x,y)$ del corpo in compressione ( $\theta_{geom} \pm \alpha > 0$)   è quello a valle di un onda d' urto di un cuneo avente un angolo $\delta_{cuneo} = atan( \frac {dy} {dx})$  posto nelle stesse condizioni asintotiche del corpo originario.  Il metodo ha validità se per ogni punto del corpo $\delta_{cuneo} \leq \delta_{lim}$  ed il bordo d' attacco è a spigolo vivo.
Per i punti del corpo in espansione ($\theta_{geom} \pm \alpha \leq 0$) si considera invece un espansione alla Prandtl e Meyer.

Nel limite della teoria dei piccoli disturbi si ha che il $c_p$ sia sul dorso che sul ventre del profilo è funzione del parametro di similitudine $K=M_\infty \theta$  (dove $\theta$ è l'angolo di deviazione della corrente, $\theta= \theta_{geom}\pm\alpha$) secondo la eq. (2)
$$
c_p= \frac {2} {M_\infty^2} \big[K + \frac{n+1}{2n}K^2 \big] + \mathcal O(K^3)  \tag 2
$$

Tuttavia per il profilo dell F-104 quest' approssimazione non è applicabile in quanto:

- sul dorso l' angolo di deviazione massimo è 25°,  in particolare tale deviazione si ha sul bordo d' uscita
- sul ventre l' angolo di deviazione massimo è sempre 25* (dato che il profilo è simmetrico), in particolare tale deviazione si ha sul bordo d' attacco

da queste considerazioni per applicare la teoria del cono tangente bisogna considerare compressioni ed espansioni non dovute ad onde di mach ma  rispettivamente dovute a compressioni da onde d' urto se l' angolo di deviazione  della corrente $\theta=\theta_{geom}+\alpha$  è positivo mentre espansioni alla prandtl e meyer se $\theta$  è negativo.
<!--TODO insert angle image-->


### Metodo urto-espansione

Nel metodo urto-espansione si tiene conto dell' interazione tra l' onda d' urto che si genera sul bordo d' attacco con le onda di espansione a valle di essa sul profilo.
In particolare passando dal profilo continuo a quello discretizzato, da variazioni della corrente infinitesime a finite:

- l' onda d'urto curva diventa un onda d' urto obliqua a tratti
- l' espansione continua alla Prandtl e Meyer è sostituita da una successione di ventagli di espansione in corrispondenza di ogni spigolo in cui c'è variazione della direzione della corrente
- I ventagli di espansione possono essere pensati come concentrati in una singola onda di espansione data la piccolezza della deviazione della corrente $\Delta\theta$

L' approssimazione fondamentale del metodo urto-espansione consiste nel **non** considerare le successive famiglie di onde riflesse nell' interazione tra onda d'urto e onda di espansione.

Ai fini del caloclo del $c_p$ sul profilo basta considerare il primo urto sul bordo d'attacco e le successive espansioni dovute alla deviazione della corrente $\delta= \theta-\theta_1$ , dove $\theta_1$ è l' angolo della corrente deviata dall' unica onda d' urto. Dunque:
$$
\begin{align}
\frac p {p_\infty}&= \frac {p_2}{p_\infty}\big|_{urto} \frac p {p_2} \big|_{espansione}\\
\frac {p_2}{p_\infty}&= \frac 1 {n+1} \big((n+2)Ma_\infty^2\sin^2(\beta) -1 \big)\\
\frac p {p_2}&= \big(1+ \big(\frac {Ma_2\delta} n\big) \big)^{n\gamma}
\end{align}
$$
Dove $\frac {p_2}{p_\infty}$ è calcolato una sola volta,  e $\beta$  è l' angolo d'urto relativo alla deviazione della corrente sul bordo d' attacco.

### Risultati



## Calcolo dei coefficienti aerodinamici

Per ogni metodo è possibile ottenere la distribuzione di $c_p$ sul profilo grazie alla quale è possibile calcolare i coefficienti di portanza $c_l$ e resistenza $c_d$  noto l' angolo d' attacco. Infatti una volta fissato il sistema di riferimento del profilo si ha che 
$$
\begin{align}

c_{F_x}&= \int_\gamma c_p \boldsymbol n\cdot \boldsymbol i dl =  \int_\gamma c_p \cos(\theta \pm \frac \pi 2) dl \\
c_{F_y}&= \int_\gamma c_p \boldsymbol n\cdot \boldsymbol j dl = \int_\gamma c_p \sin(\theta \pm \frac \pi 2) dl  \\

\end{align}
$$
da cui si possono calcolare i coefficienti di portanza e resistenza

$$
\begin{align}
c_d&= c_{F_x}\cos(\alpha) + c_{F_y}\sin(\alpha) \\
c_l&= -c_{F_x}\sin(\alpha) + c_{F_y}\cos(\alpha)\\
\end{align}
$$

nel codice il calcolo dei coefficienti aerodinamici è implementato attraverso la routine <code>CalcClCd()</code> , gli integrali per il calcolo di $C_{F_x}$ e $C_{F_y}$ sono approssimati attraverso una semplice formula di quadratura.  Per esempio per il calcolo di $C_{F_x}$ si ha:
$$
C_{F_x}= \sum _{dorso} c_{p_i} \Delta l_i \cos(\theta + \frac \pi 2) - \sum_{ventre} c_{p_i} \Delta l_i\cos(\theta + \frac \pi 2)
$$

