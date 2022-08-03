import { Component, Inject, ElementRef, OnInit } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
import * as d3 from 'd3';

export const digitRound = (x: number, ff = 1e5) => {
  const xff = x * ff;
  const cxff = Math.ceil(xff);
  let delt = xff - cxff;
  if (Math.abs(delt) < 1e-8) {
    return cxff / ff;
  }
  const fxff = Math.floor(xff);
  delt = xff - fxff;
  if (Math.abs(delt) < 1e-8) {
    return fxff / ff;
  }
  return x;
}
@Component({
  selector: 'app-optimisegeneral',
  templateUrl: './optimisegeneral.component.html',
  styleUrls: ['./optimisegeneral.component.css']
})

export class OptimisegeneralComponent implements OnInit {
  width = 1100;//width and height for weight graph
  height = 300;
  generalfile = "GLdist";
  filterzero = false;
  shortside = false;
  format = d3.format('0.6f')
  opt: Optimise = {} as Optimise;
  constructor(private http: HttpClient, @Inject('BASE_URL') private baseUrl: string, public element: ElementRef) {
    http.get<Optimise>(baseUrl + 'optimise/general?doOpt=false&datafile=' + this.generalfile).subscribe(result => {
      this.opt = result;
      d3.select(this.element.nativeElement).select('input.checkzero').attr('value', this.filterzero);
      console.log(this.opt, result);
      this.shortside = false;
      for (let i = 0; i < this.opt.n; ++i) {
        if (this.opt.l[i] < 0) this.shortside = true;
      }
      console.log(this.shortside);
    }, error => console.error(error));
  }
  over(e: MouseEvent, inout = false) {
    d3.select(e.target as HTMLInputElement & EventTarget).classed('over', inout);
    //console.log(e.target);
  }
  ngOnInit(): void {
    d3.select(this.element.nativeElement).select('input.checkzero').attr('value', this.filterzero);
    d3.select(this.element.nativeElement).select('input.filer').attr('value', this.generalfile);
  }
  sendStep() {
    this.opt.doOpt = true;
    let back: string | number | null | undefined = +(d3.select(this.element.nativeElement).select('input.step.g').node() as HTMLInputElement & Event).value;
    this.opt.gamma = back;
    try { back = +(d3.select(this.element.nativeElement).select('input.step.round').node() as HTMLInputElement & Event).value; } catch { back = null; }
    if (back != undefined) this.opt.round = back;
    back = +(d3.select(this.element.nativeElement).select('input.step.delt').node() as HTMLInputElement & Event).value;
    this.opt.delta = back;
    back = +(d3.select(this.element.nativeElement).select('input.step.basket').node() as HTMLInputElement & Event).value;
    this.opt.basket = back;
    back = +(d3.select(this.element.nativeElement).select('input.step.trades').node() as HTMLInputElement & Event).value;
    this.opt.trades = back;
    back = +(d3.select(this.element.nativeElement).select('input.step.minhold').node() as HTMLInputElement & Event).value;
    this.opt.min_holding = back;
    back = +(d3.select(this.element.nativeElement).select('input.step.mintrade').node() as HTMLInputElement & Event).value;
    this.opt.min_trade = back;
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
    try { back = +(d3.select(this.element.nativeElement).select('input.step.VU').node() as HTMLInputElement & Event).value; } catch { back = null; }
    if (back != undefined) this.opt.value = back;
    if (back != undefined) try { back = +(d3.select(this.element.nativeElement).select('input.step.VL').node() as HTMLInputElement & Event).value; } catch { back = null; }
    if (back != undefined) this.opt.valuel = back;
    try { back = +(d3.select(this.element.nativeElement).select('input.step.RMU').node() as HTMLInputElement & Event).value; } catch { back = null; }
    if (back != undefined) this.opt.rmax = back;
    if (back != undefined) try { back = +(d3.select(this.element.nativeElement).select('input.step.RML').node() as HTMLInputElement & Event).value; } catch { back = null; }
    if (back != undefined) this.opt.rmin = back;
    if (this.opt.losSmax != null && this.opt.losSmin != null) this.opt.losSopt = true;
    else if (this.opt.etLmax != null && this.opt.etLmin != null) this.opt.etLopt = true;
    this.opt.logfile = "bitapluslog";
    this.sendData('optimise/general', this.opt)
      .subscribe(ddd => {
        console.log(ddd);
        this.opt = ddd;
        this.opt.w = this.opt.w.map(x => digitRound(x));
        console.log(this.opt.w);
        this.shortside = false;
        d3.select(this.element.nativeElement).select('input.checkzero').attr('value', this.filterzero);
        for (let i = 0; i < this.opt.n; ++i) {
          if (this.opt.l[i] < 0) this.shortside = true;
        }
        console.log(this.shortside);
      }, error => console.error(error));
  }
  zerooff(e: Event) {
    this.filterzero = !this.filterzero;
    d3.select(this.element.nativeElement).select('input.checkzero').attr('value', this.filterzero);
    console.log((d3.select(e.target as HTMLInputElement & EventTarget).node() as HTMLInputElement & Event).value);
    console.log((d3.select(this.element.nativeElement).select('input.checkzero').attr('value')));
    //   this.opt=this.opt;
  }
  file(e: Event) {
    console.log(e.target);//Must use .node() to get the updated value
    //   console.log((d3.select(this.element.nativeElement).select('input.filer').node() as HTMLInputElement & Event).value);
    this.generalfile = (d3.select(e.target as HTMLInputElement & EventTarget).node() as HTMLInputElement & Event).value;
    //console.log(d3.select(e.target as HTMLInputElement & EventTarget).attr('value'));
    //console.log(d3.select(e.target as HTMLInputElement).attr('value'));
    //   console.log(this.generalfile);
    this.http.get<Optimise>(this.baseUrl + 'optimise/general?doOpt=false&datafile=' + this.generalfile).subscribe(result => {
      this.opt = result;
      d3.select(this.element.nativeElement).select('input.checkzero').attr('value', this.filterzero);
      console.log(this.opt, result);
      this.shortside = false;
      for (let i = 0; i < this.opt.n; ++i) {
        if (this.opt.l[i] < 0) this.shortside = true;
      }
      console.log(this.shortside);
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
  cval:Array<number>,
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
  versionString:string,
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
  datafile: string,
  logfile: string
}