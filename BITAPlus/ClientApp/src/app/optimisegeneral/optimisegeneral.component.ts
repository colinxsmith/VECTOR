import { Component, Inject, ElementRef } from '@angular/core';
import { HttpClient, HttpHeaders, HttpParams } from '@angular/common/http';
import * as d3 from 'd3';

@Component({
  selector: 'app-optimisegeneral',
  templateUrl: './optimisegeneral.component.html',
  styleUrls: ['./optimisegeneral.component.css']
})

export class OptimisegeneralComponent {
  width = 1000;//width and height for weight graph
  height = 300;
  format = d3.format('0.6f')
  opt: Optimise = {} as Optimise;
  constructor(private http: HttpClient, @Inject('BASE_URL') private baseUrl: string, public element: ElementRef) {
    http.get<Optimise>(baseUrl + 'optimise/general').subscribe(result => {
      this.opt = result;
      console.log(this.opt, result);
    }, error => console.error(error));
  }
  over(e: MouseEvent, inout = false) {
    d3.select(e.target as HTMLInputElement & EventTarget).classed('over', inout);
  }
  sendStep() {
    let back: string | number | null | undefined = +(d3.select(this.element.nativeElement).select('input.step.g').node() as HTMLInputElement & Event).value;
    this.opt.gamma = back;
    back = +(d3.select(this.element.nativeElement).select('input.step.delt').node() as HTMLInputElement & Event).value;
    this.opt.delta = back;
    try { back = +(d3.select(this.element.nativeElement).select('input.step.EU').node() as HTMLInputElement & Event).value; } catch { back = null; }
    if (back != undefined) this.opt.etLmax = back;
    try { back = +(d3.select(this.element.nativeElement).select('input.step.EL').node() as HTMLInputElement & Event).value; } catch { back = null; }
    if (back != undefined) this.opt.etLmin = back;
    try { back = +(d3.select(this.element.nativeElement).select('input.step.LU').node() as HTMLInputElement & Event).value; } catch { back = null; }
    if (back != undefined) this.opt.losSmax = back;
    try { back = +(d3.select(this.element.nativeElement).select('input.step.LL').node() as HTMLInputElement & Event).value; } catch { back = null; }
    if (back != undefined) this.opt.losSmin = back;
    try { back = +(d3.select(this.element.nativeElement).select('input.step.RU').node() as HTMLInputElement & Event).value; } catch { back = null; }
    if (back != undefined) this.opt.maxRisk = back;
    if (back != undefined) try { back = +(d3.select(this.element.nativeElement).select('input.step.RL').node() as HTMLInputElement & Event).value; } catch { back = null; }
    if (back != undefined) this.opt.minRisk = back;
    if (this.opt.losSmax != null && this.opt.losSmin != null) this.opt.losSopt = true;
    else if (this.opt.etLmax != null && this.opt.etLmin != null) this.opt.etLopt = true;
    this.sendData('optimise/general', this.opt)
      .subscribe(ddd => {
        console.log(ddd);
        this.opt = ddd;
        console.log(this.opt.w);
      }, error => console.error(error));
  }
  sendData(key = 'optimise/general', sendObject: Optimise) {
    const options = {
      headers: new HttpHeaders()
        .set('Content-Type', 'application/json')
    };
    return this.http.post<Optimise>(`${this.baseUrl}${key}`, sendObject, options);
  }
}

interface Result {


  gross: number,
  longvalue: number,
  shortvalue: number,
  shortoverlong: number,
  minhold: number,
  mintrade: number,
  cost: number,
  turnover: number,
  basket: number,
  trades: number,
  var: number,
  vaRindex: number,
  etl: number,
  loss: number,
  breakdown: Array<number>,
  risk: number,
  expreturn: number,
  mctr: Array<number>,
  fmctr: Array<number>,
  sPmctr: Array<number>,
  fx: Array<number>,
  facrisk: number,
  specrisk: number
}
interface Optimise {
  doOpt: boolean,
  result: Result,
  back: number,
  message: string,
  ogamma: number,
  minRisk: number,
  maxRisk: number,
  rmax: number,
  rmin: number,
  min_holding: number,
  min_trade: number,
  basket: number,
  trades: number,
  tlen: number,
  nfac: number,
  nabs: number,
  mabs: number,
  n: number,
  m: number,
  bench: Array<number>,
  l: Array<number>,
  u: Array<number>,
  a: Array<number>,
  data: Array<number>,
  buy: Array<number>,
  abs_A: Array<number>,
  abs_L: Array<number>,
  abs_U: Array<number>,
  mask: Array<number>,
  longbasket: number,
  shortbasket: number,
  tradesell: number,
  tradebuy: number,
  sv: Array<number>,
  fc: Array<number>,
  fl: Array<number>,
  i_A: Array<number>,
  sell: Array<number>,
  delta: number,
  gamma: number,
  value: number,
  valuel: number,
  kappa: number,
  gstrength: number,
  tail: number,
  etLopt: boolean,
  etLmin: number,
  etLmax: number,
  names: Array<string>,
  initial: Array<number>,
  q: Array<number>,
  w: Array<number>,
  min_lot: Array<number>,
  size_lot: Array<number>,
  shake: Array<number>,
  alpha: Array<number>,
  cvargLprob: boolean,
  targetReturn: Array<number>,
  losSopt: boolean,
  losSmin: number,
  losSmax: number,
  round: number,
  datafile: string
}